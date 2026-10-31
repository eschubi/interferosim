# src/interferosim/signal/iq.py -- simulated SDR receiver signals (interferosim.signal.iq)
"""Simulation of the IQ samples the SDRs of an interferometer record.

Signal model
------------
Every sky component radiates noise: an unresolved source is a stationary,
complex Gaussian signal s(t) whose power is its flux density S. A station
k at position r_k receives it delayed by the geometric delay
tau_k(t) = -((r_k - r_ref) . s) / c relative to the reference station and,
after mixing to complex baseband with the LO at f_c, rotated by the phase
-2 pi f_c tau_k(t):

    x_k(t) = sum_sources sqrt(P_k(theta)) s(t - tau_k) exp(-2 pi i f_c tau_k)
             + n_k(t),

P_k is the primary beam of the antenna of station k at the angle theta of
the source from the pointing direction, n_k the receiver noise of power
SEFD_k (independent per station). Units: sqrt(Jy) per sample, so that
``<|x_k|^2> = SEFD_k + sum S P_k`` and, for two stations,
``<x_i x_j*> = sum S sqrt(P_i P_j) exp(-2 pi i (w + u l + v m))``.

Extended components (Gaussian, disc) are incoherent sums of such signals.
The simulator does not split them into points: per block it builds the
covariance matrix of the station signals directly from the analytic
component visibilities, R_ij = sum_c S_c sqrt(P_i P_j) V_c(u_ij, v_ij) / S_c
(R_kk = sum_c S_c P_k), draws correlated white streams with this covariance
and then applies the delay and fringe phase of the phase centre to every
station. The expectation of the simulated correlation is therefore exactly
the visibility of the sky model (``expected_visibility``).

Instrument errors per station (:class:`StationErrors`): LO frequency offset
and phase (frequency reference), clock offset (time reference, e.g. PPS);
they turn into a residual fringe rotation and a residual delay in the
correlation. :func:`quantize` imitates an ADC / host format of a few bits.

Origin
------
AAS6 ``modules/sdr/iq_record_node.py`` (``_SimDualSource``: common sky
signal, station noise) and ``modules/interfero/skysim.py`` (delay and
fringe phase of the non-reference station, primary beam); extended here to
several sources, any number of stations and instrument errors.

Version 0.4.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np

from ..constants import C_LIGHT
from ..instrument import Array, Receiver
from ..observation.geometry import (SIDEREAL_DEG_PER_S, direction, elevation_deg,
                                    geometric_delay, gha_deg, uvw)
from ..sky import SkyModel
from .correlate import fractional_delay

__all__ = ["IQBlock", "StationErrors", "IQSimulator", "quantize"]


# ======================================================================
@dataclass
class IQBlock:
    """Complex baseband samples of one station."""

    station: str
    data: np.ndarray                 # complex64, shape (n,)
    sample_rate_sps: float
    center_freq_hz: float
    t0_s: float = 0.0                # start time relative to the simulation start
    meta: Dict[str, object] = field(default_factory=dict)

    @property
    def duration_s(self) -> float:
        return self.data.size / self.sample_rate_sps

    def times(self) -> np.ndarray:
        """Sample times in seconds from the simulation start."""
        return self.t0_s + np.arange(self.data.size) / self.sample_rate_sps

    def save_sigmf(self, path: Union[str, Path]) -> Path:
        """Write the block as SigMF recording: ``<path>.sigmf-data``
        (complex float32, little endian) and ``<path>.sigmf-meta`` (JSON).
        Returns the path of the meta file."""
        base = Path(path)
        base.parent.mkdir(parents=True, exist_ok=True)
        data_path = base.with_suffix(".sigmf-data")
        self.data.astype("<c8").tofile(data_path)
        meta = {
            "global": {"core:datatype": "cf32_le", "core:sample_rate": self.sample_rate_sps,
                       "core:version": "1.0.0", "core:recorder": "interferosim",
                       "core:description": f"simulated IQ samples, station {self.station}",
                       "interferosim:meta": {k: v for k, v in self.meta.items()
                                             if isinstance(v, (int, float, str, bool, list))}},
            "captures": [{"core:sample_start": 0, "core:frequency": self.center_freq_hz,
                          "interferosim:t0_s": self.t0_s}],
            "annotations": [],
        }
        meta_path = base.with_suffix(".sigmf-meta")
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return meta_path


# ======================================================================
@dataclass
class StationErrors:
    """Instrument errors of one station (all zero: ideal station)."""

    lo_offset_hz: float = 0.0        # frequency error of the LO
    lo_phase_rad: float = 0.0        # phase of the LO at t = 0
    clock_offset_s: float = 0.0      # sample clock / time stamp error

    @classmethod
    def from_receiver(cls, receiver: Receiver, rng: np.random.Generator,
                      clock_source: Optional[str] = None,
                      f_hz: Optional[float] = None) -> "StationErrors":
        """Random errors within the data-sheet limits of ``receiver``:
        LO offset uniform in +-y f, LO phase uniform, clock offset uniform in
        +-PPS accuracy if the time comes from a GPSDO (else 0)."""
        df = receiver.lo_offset_hz(f_hz, clock_source) or 0.0
        dt = 0.0
        if receiver.time_source == "gpsdo" and receiver.pps_accuracy_ns:
            dt = receiver.pps_accuracy_ns * 1e-9
        return cls(lo_offset_hz=float(rng.uniform(-df, df)),
                   lo_phase_rad=float(rng.uniform(0, 2 * math.pi)),
                   clock_offset_s=float(rng.uniform(-dt, dt)))


# ======================================================================
def quantize(x: np.ndarray, bits: int, rms_fraction: Optional[float] = None) -> np.ndarray:
    """Quantise I and Q to ``bits`` bits (two's complement levels) and return
    the result in the original units (levels times step).

    Without ``rms_fraction`` the step is the optimum for Gaussian noise
    (Max 1960: 0.996, 0.586, 0.335 sigma for 2, 3, 4 bits); from 5 bits on the
    signal rms is put at 1/4 of full scale. With ``rms_fraction`` the rms is
    put at that fraction of full scale. 1 bit keeps the sign only."""
    if bits < 1:
        raise ValueError("bits must be >= 1")
    rms = float(np.sqrt(np.mean(np.abs(x) ** 2) / 2.0)) or 1.0
    if bits == 1:                                      # sign only
        return (np.sign(x.real) + 1j * np.sign(x.imag)).astype(np.complex64) * rms
    full = 2 ** (bits - 1)
    if rms_fraction is None:
        step = {2: 0.996, 3: 0.586, 4: 0.335}.get(bits, 1.0 / (0.25 * full)) * rms
    else:
        step = rms / (rms_fraction * full)

    def q(v):
        # mid-rise quantiser: levels (k + 1/2) step, k = -full .. full-1
        return (np.clip(np.floor(v / step), -full, full - 1) + 0.5) * step
    return (q(x.real) + 1j * q(x.imag)).astype(np.complex64)


# ======================================================================
class IQSimulator:
    """Simulated IQ recording of all stations of an array.

    Parameters
    ----------
    sky
        Sky model; its field centre (RA/Dec) is the pointing and phase centre.
    array
        Stations (positions, antennas); the first station is the reference.
    receiver
        Receiver (sample rate, centre frequency, T_sys); default
        ``array.receiver``. ``sample_rate_sps`` / ``center_freq_hz`` override.
    errors
        Instrument errors per station name (default: ideal).
    seed
        Seed of the random generators (sky signal, noise, sub-sources).
    """

    def __init__(self, sky: SkyModel, array: Array, receiver: Optional[Receiver] = None, *,
                 sample_rate_sps: Optional[float] = None,
                 center_freq_hz: Optional[float] = None,
                 errors: Optional[Dict[str, StationErrors]] = None,
                 seed: int = 1, block_size: int = 16384,
                 receiver_noise: bool = True, sky_signal: bool = True):
        if sky.center.ra_deg is None or sky.center.dec_deg is None:
            raise ValueError("the sky model needs a field centre with RA/Dec")
        self.sky = sky
        self.array = array
        self.receiver = receiver or array.receiver
        if self.receiver is None:
            raise ValueError("no receiver: pass one or use an array with a receiver")
        self.fs = float(sample_rate_sps or self.receiver.sample_rate_sps)
        self.fc = float(center_freq_hz or self.receiver.center_freq_hz or 1.420405751768e9)
        self.errors = errors or {}
        self.seed = int(seed)
        self.block_size = int(block_size)
        self.receiver_noise = receiver_noise
        self.sky_signal = sky_signal

    # ------------------------------------------------------------------
    @property
    def wavelength_m(self) -> float:
        return C_LIGHT / self.fc

    # ------------------------------------------------------------------
    def gha_deg(self, t_s, lha0_deg: float = 0.0) -> np.ndarray:
        """Greenwich hour angle at time ``t_s`` after the start, for the
        local hour angle ``lha0_deg`` at the reference station at t = 0."""
        ref = self.array.reference
        return gha_deg(self.sky.center.ra_deg, lha_deg=lha0_deg + SIDEREAL_DEG_PER_S *
                       np.asarray(t_s, dtype=float), lon_deg=ref.lon_deg)

    def uvw(self, t_s, lha0_deg: float = 0.0) -> np.ndarray:
        """(u, v, w) of all baselines, shape (n_baselines, n_t, 3)."""
        g = self.gha_deg(t_s, lha0_deg)
        return np.stack([uvw(b, self.sky.center.dec_deg, g, self.wavelength_m)
                         for b in self.array.baselines()])

    def phase_centre_delay_s(self, station: int, t_s, lha0_deg: float = 0.0) -> np.ndarray:
        """Geometric delay of the phase centre for station ``station``
        relative to the reference, at times ``t_s``."""
        r = [s.ecef_m() for s in self.array.stations]
        g = self.gha_deg(t_s, lha0_deg)
        return geometric_delay(r[station] - r[0], direction(self.sky.center.dec_deg, g))

    def elevation_deg(self, t_s, lha0_deg: float = 0.0) -> np.ndarray:
        """Elevation of the phase centre at the reference station."""
        ref = self.array.reference
        g = self.gha_deg(t_s, lha0_deg)
        return elevation_deg(direction(self.sky.center.dec_deg, g), ref.lat_deg, ref.lon_deg)

    def _beam_amplitudes(self) -> np.ndarray:
        """sqrt(P_k) of every component at its centre, shape (n_st, n_comp)."""
        out = np.ones((len(self.array.stations), len(self.sky)))
        for k, st in enumerate(self.array.stations):
            if st.antenna is None:
                continue
            for c_i, c in enumerate(self.sky):
                th = math.hypot(c.dx_deg, c.dy_deg)
                out[k, c_i] = math.sqrt(max(float(st.antenna.beam(th, self.fc)), 0.0))
        return out

    def expected_visibility(self, i: int, j: int, u, v) -> np.ndarray:
        """Expectation of <x_i x_j*> after delay correction and fringe
        stopping: sum_c sqrt(P_i P_j) V_c(u, v) (primary beam at the centre of
        each component)."""
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        a = self._beam_amplitudes()
        V = np.zeros(np.broadcast(u, v).shape, complex)
        for c_i, c in enumerate(self.sky):
            V = V + a[i, c_i] * a[j, c_i] * c.visibility(u, v)
        return V

    def station_sky_power_jy(self, k: int) -> float:
        """Sky power received by station k, sum_c S_c P_k."""
        a = self._beam_amplitudes()
        return float(sum(c.flux_jy * a[k, c_i] ** 2 for c_i, c in enumerate(self.sky)))

    def sigma_visibility_jy(self, i: int, j: int, t_avg_s: float) -> float:
        """Expected noise of the real (imaginary) part of one averaged
        visibility, sqrt((SEFD_i + S_i)(SEFD_j + S_j) / (2 N)), N = fs t_avg."""
        st = self.array.stations
        tot = []
        for k in (i, j):
            sefd = self.receiver.sefd_jy(st[k].antenna, self.fc) if st[k].antenna else 0.0
            tot.append(sefd * self.receiver_noise + self.station_sky_power_jy(k) * self.sky_signal)
        return math.sqrt(tot[0] * tot[1] / (2.0 * self.fs * t_avg_s))

    # ------------------------------------------------------------------
    def simulate(self, duration_s: float, lha0_deg: float = 0.0) -> Dict[str, IQBlock]:
        """IQ samples of all stations for ``duration_s`` seconds, starting at
        the local hour angle ``lha0_deg`` of the phase centre at the
        reference station. Returns {station name: IQBlock}.

        Method (per block of ``block_size`` samples): the sky signal of all
        point sources is generated as n_station correlated white streams with
        the covariance R = C C^H, C[k, j] = sqrt(S_j P_k) exp(-2 pi i f_c
        dtau_kj), where dtau_kj is the delay of source j relative to the
        phase centre; each stream is then delayed by the phase-centre delay
        tau0_k (plus clock offset) and rotated by exp(-2 pi i f_c tau0_k(t))
        sample by sample. This is exact as long as the differential delays
        dtau are small against one sample (no bandwidth smearing; checked,
        see ``meta['max_differential_delay_samples']``).
        """
        n_total = int(round(duration_s * self.fs))
        if n_total <= 0:
            raise ValueError("duration too short")
        stations = self.array.stations
        n_st = len(stations)
        r = np.array([s.ecef_m() for s in stations])
        bvec = r - r[0]
        dec = self.sky.center.dec_deg
        rec = self.receiver
        sefd = [rec.sefd_jy(s.antenna, self.fc) if s.antenna else 0.0 for s in stations]
        amp = self._beam_amplitudes()                               # (n_st, n_comp)
        pairs = [(i, j) for i in range(n_st) for j in range(i + 1, n_st)]
        # largest differential delay (component edge) against one sample
        reach = max((math.radians(math.hypot(c.dx_deg, c.dy_deg) + 0.5 * c.extent_deg)
                     for c in self.sky), default=0.0)
        max_dd = float(np.linalg.norm(bvec, axis=1).max() * reach / C_LIGHT * self.fs)
        err = [self.errors.get(s.name, StationErrors()) for s in stations]

        rng_sig = np.random.default_rng(self.seed + 1)
        rng_noise = [np.random.default_rng(self.seed + 100 + k) for k in range(n_st)]
        out = [np.empty(n_total, np.complex64) for _ in stations]
        max_delay = np.linalg.norm(bvec, axis=1).max() / C_LIGHT
        max_clock = max(abs(e.clock_offset_s) for e in err)
        pad = int(math.ceil((max_delay + max_clock) * self.fs)) + 8
        g_eps = 0.01                                                # degrees, for the delay rate

        for b0 in range(0, n_total, self.block_size):
            n = min(self.block_size, n_total - b0)
            t = (b0 + np.arange(n)) / self.fs
            t_mid = t[0] + 0.5 * n / self.fs
            g_mid = float(self.gha_deg(t_mid, lha0_deg))
            s0 = direction(dec, g_mid)
            tau0 = -(bvec @ s0) / C_LIGHT                                   # (n_st,)
            tau0b = -(bvec @ direction(dec, g_mid + g_eps)) / C_LIGHT
            rate0 = (tau0b - tau0) / (g_eps / SIDEREAL_DEG_PER_S)
            if self.sky_signal and len(self.sky):
                # covariance of the station sky signals from the component visibilities
                R = np.zeros((n_st, n_st), complex)
                for k in range(n_st):
                    R[k, k] = sum(c.flux_jy * amp[k, c_i] ** 2 for c_i, c in enumerate(self.sky))
                for (i, j) in pairs:
                    u_, v_, _ = uvw(bvec[j] - bvec[i], dec, g_mid, self.wavelength_m)
                    vij = sum(amp[i, c_i] * amp[j, c_i] * c.visibility(u_, v_)
                              for c_i, c in enumerate(self.sky))
                    R[i, j] = complex(vij)
                    R[j, i] = np.conj(R[i, j])
                w, U = np.linalg.eigh(R)
                L = U * np.sqrt(np.clip(w, 0.0, None))[None, :]
                Z = (rng_sig.standard_normal((n_st, n + 2 * pad))
                     + 1j * rng_sig.standard_normal((n_st, n + 2 * pad))) / math.sqrt(2.0)
                Y = L @ Z
            else:
                Y = np.zeros((n_st, n + 2 * pad), complex)
            for k, e in enumerate(err):
                d = (tau0[k] + e.clock_offset_s) * self.fs
                y = fractional_delay(Y[k], d) if d else Y[k]
                y = y[pad:pad + n]
                tau_t = tau0[k] + rate0[k] * (t - t_mid)
                if k > 0 or tau0[k]:
                    y = y * np.exp(-2j * np.pi * self.fc * tau_t)
                if self.receiver_noise and sefd[k] > 0:
                    y = y + (rng_noise[k].standard_normal(n) + 1j * rng_noise[k].standard_normal(n)) \
                        * math.sqrt(sefd[k] / 2.0)
                if e.lo_offset_hz or e.lo_phase_rad:
                    y = y * np.exp(1j * (2 * np.pi * e.lo_offset_hz * t + e.lo_phase_rad))
                out[k][b0:b0 + n] = y.astype(np.complex64)

        el = self.elevation_deg(np.array([0.0, duration_s]), lha0_deg)
        blocks = {}
        for k, s in enumerate(stations):
            blocks[s.name] = IQBlock(
                station=s.name, data=out[k], sample_rate_sps=self.fs, center_freq_hz=self.fc,
                meta={"sky": self.sky.name, "ra_deg": self.sky.center.ra_deg, "dec_deg": dec,
                      "lha0_deg": lha0_deg, "sefd_jy": sefd[k],
                      "elevation_start_deg": float(el[0]), "elevation_end_deg": float(el[1]),
                      "lo_offset_hz": err[k].lo_offset_hz, "clock_offset_s": err[k].clock_offset_s,
                      "seed": self.seed,
                      "max_differential_delay_samples": max_dd})
        return blocks
