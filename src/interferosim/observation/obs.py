# src/interferosim/observation/obs.py -- observation, uv-coverage, visibility and transit simulation (interferosim.observation.obs)
"""Earth rotation synthesis: observations, uv-coverage, simulated visibilities.

* :class:`Observation` -- an array observing a direction (RA/Dec) over a
  range of hour angles on a given date: hour angles, times, elevations per
  station, uvw per baseline, geometric delay and fringe rate, the mask of
  usable integrations (all stations above ``min_elev_deg``).
* :class:`UVCoverage` -- the uv-tracks of an observation.
* :class:`VisibilitySimulator` -- analytic visibilities of a sky model with
  radiometer noise (tracking antennas, fringe-stopped), as
  :class:`~interferosim.visdata.VisData`.
* :class:`TransitSimulator` -- drift scan of fixed antennas (notebook B,
  Sun with two horns): the source moves through the primary beams, the
  phase advances with the fringe rate (no fringe stopping), noise per dump.

Origin: ``aas6sim`` 1.3.1 (``model.uv_track_ew``, ``model.elevation_deg``,
``simulate.solar_transit``); here for any baseline orientation with the
ECEF geometry of :mod:`interferosim.observation.geometry`.

Version 0.5.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from ..constants import C_LIGHT
from ..instrument import Array, Receiver
from ..sky import SkyModel
from ..visdata import VisData
from .geometry import (SIDEREAL_DEG_PER_S, delay_rate, direction, elevation_deg,
                       geometric_delay, gmst_deg, sky_basis, uvw)
from .sun import mjd_from_iso

__all__ = ["Observation", "UVCoverage", "VisibilitySimulator", "TransitSimulator"]

SOLAR_PER_SIDEREAL = 0.99726957          # length of a sidereal day in solar days


# ======================================================================
@dataclass
class Observation:
    """An array observing the direction (ra_deg, dec_deg)."""

    array: Array
    ra_deg: float
    dec_deg: float
    freq_hz: float = 1.420405751768e9
    ha_range_h: Tuple[float, float] = (-6.0, 6.0)   # hour angle at the reference station
    t_int_s: float = 60.0                            # integration (dump) time
    min_elev_deg: float = 10.0
    date: str = "2026-11-01"                          # UTC date of the transit used
    name: str = ""

    @classmethod
    def from_sky(cls, sky: SkyModel, array: Array, **kw) -> "Observation":
        if sky.center.ra_deg is None or sky.center.dec_deg is None:
            raise ValueError("the sky model needs a field centre with RA/Dec")
        return cls(array=array, ra_deg=sky.center.ra_deg, dec_deg=sky.center.dec_deg,
                   name=kw.pop("name", sky.center.label or sky.name), **kw)

    # ------------------------------------------------------------------
    @property
    def wavelength_m(self) -> float:
        return C_LIGHT / self.freq_hz

    def hour_angles(self) -> np.ndarray:
        """Hour angles [deg] at the reference station, one per integration
        (centre of each dump)."""
        h0, h1 = self.ha_range_h
        dh = SIDEREAL_DEG_PER_S * self.t_int_s
        n = max(int(round((h1 - h0) * 15.0 / dh)), 1)
        return h0 * 15.0 + (np.arange(n) + 0.5) * dh

    def transit_mjd(self) -> float:
        """MJD of the upper transit at the reference station on ``date``."""
        mjd0 = mjd_from_iso(self.date)
        lon = self.array.reference.lon_deg
        lst0 = float(gmst_deg((mjd0 - 40587.0) * 86400.0)) + lon
        dh = (self.ra_deg - lst0) % 360.0                 # sidereal degrees to go
        return mjd0 + dh / 360.0 * SOLAR_PER_SIDEREAL

    def times_mjd(self) -> np.ndarray:
        return self.transit_mjd() + self.hour_angles() / 360.0 * SOLAR_PER_SIDEREAL

    def gha_deg(self, ha_deg=None) -> np.ndarray:
        h = self.hour_angles() if ha_deg is None else np.asarray(ha_deg, dtype=float)
        return np.mod(h - self.array.reference.lon_deg, 360.0)

    def elevation_deg(self, ha_deg=None) -> np.ndarray:
        """Elevation of the source at every station, shape (n_st, n_t)."""
        s = direction(self.dec_deg, self.gha_deg(ha_deg))
        return np.array([elevation_deg(s, st.lat_deg, st.lon_deg) for st in self.array.stations])

    def visible(self, ha_deg=None) -> np.ndarray:
        """True where the source is above ``min_elev_deg`` at all stations."""
        return np.all(self.elevation_deg(ha_deg) >= self.min_elev_deg, axis=0)

    def uvw(self, ha_deg=None) -> np.ndarray:
        """(u, v, w) of all baselines, shape (n_bl, n_t, 3), in wavelengths."""
        g = self.gha_deg(ha_deg)
        return np.stack([uvw(b, self.dec_deg, g, self.wavelength_m) for b in self.array.baselines()])

    def delay_s(self, ha_deg=None) -> np.ndarray:
        """Geometric delay of the phase centre, shape (n_bl, n_t)."""
        s = direction(self.dec_deg, self.gha_deg(ha_deg))
        return np.stack([geometric_delay(b, s) for b in self.array.baselines()])

    def fringe_rate_hz(self, ha_deg=None) -> np.ndarray:
        """Natural fringe rate f_c d tau / dt, shape (n_bl, n_t)."""
        g = self.gha_deg(ha_deg)
        return np.stack([self.freq_hz * delay_rate(b, self.dec_deg, g) for b in self.array.baselines()])


# ======================================================================
@dataclass
class UVCoverage:
    """uv-tracks of an observation (wavelengths), shape (n_bl, n_t)."""

    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    ha_deg: np.ndarray
    mask: np.ndarray
    baseline_names: List[str]
    wavelength_m: float

    @classmethod
    def from_observation(cls, obs: Observation) -> "UVCoverage":
        ha = obs.hour_angles()
        q = obs.uvw(ha)
        return cls(u=q[..., 0], v=q[..., 1], w=q[..., 2], ha_deg=ha, mask=obs.visible(ha),
                   baseline_names=obs.array.baseline_names(), wavelength_m=obs.wavelength_m)

    @property
    def max_uv(self) -> float:
        m = self.mask[None, :]
        return float(np.max(np.hypot(self.u, self.v) * m)) if self.mask.any() else 0.0

    def resolution_arcsec(self) -> float:
        """Approximate angular resolution 1 / (2 q_max) [arcsec]."""
        q = self.max_uv
        return math.degrees(1.0 / (2.0 * q)) * 3600.0 if q > 0 else float("inf")


# ======================================================================
def _beam_amp(antenna, theta_deg, f_hz):
    if antenna is None:
        return np.ones_like(np.asarray(theta_deg, dtype=float))
    return np.sqrt(np.clip(antenna.beam(np.asarray(theta_deg, dtype=float), f_hz), 0.0, None))


class VisibilitySimulator:
    """Analytic visibilities of a sky model on an observation, with noise.

    Antennas track the phase centre; the primary beam weights each
    component at its offset (constant during tracking). The visibility of
    baseline i-j after fringe stopping is sum_c sqrt(P_i P_j) V_c(u, v); the
    noise of the real (imaginary) part per integration is
    sqrt((SEFD_i + S_i)(SEFD_j + S_j) / (2 B t_int)).
    """

    def __init__(self, receiver: Optional[Receiver] = None, *, add_noise: bool = True,
                 bandwidth_hz: Optional[float] = None, seed: int = 1):
        self.receiver = receiver
        self.add_noise = add_noise
        self.bandwidth_hz = bandwidth_hz
        self.seed = seed

    def run(self, sky: SkyModel, obs: Observation, only_visible: bool = True) -> VisData:
        rec = self.receiver or obs.array.receiver
        if rec is None:
            raise ValueError("no receiver")
        bw = self.bandwidth_hz or rec.bandwidth_hz
        st = obs.array.stations
        f = obs.freq_hz
        ha = obs.hour_angles()
        keep = obs.visible(ha) if only_visible else np.ones(ha.size, bool)
        ha = ha[keep]
        t_mjd = obs.times_mjd()[keep]
        q = obs.uvw(ha)
        theta = [math.hypot(c.dx_deg, c.dy_deg) for c in sky]
        amp = np.array([[float(_beam_amp(s.antenna, th, f)) for th in theta] for s in st])
        s_in = [float(sum(c.flux_jy * amp[k, ci] ** 2 for ci, c in enumerate(sky))) for k in range(len(st))]
        sefd = [rec.sefd_jy(s.antenna, f) if s.antenna else 0.0 for s in st]
        rng = np.random.default_rng(self.seed)
        rows = {k: [] for k in ("t", "ha", "u", "v", "w", "b", "vis", "wt")}
        for b_idx, (i, j) in enumerate(obs.array.baseline_pairs()):
            u, v, w = q[b_idx, :, 0], q[b_idx, :, 1], q[b_idx, :, 2]
            V = np.zeros(ha.size, complex)
            for ci, c in enumerate(sky):
                V += amp[i, ci] * amp[j, ci] * c.visibility(u, v)
            sig = math.sqrt((sefd[i] + s_in[i]) * (sefd[j] + s_in[j]) / (2.0 * bw * obs.t_int_s))
            if self.add_noise:
                V = V + sig * (rng.standard_normal(ha.size) + 1j * rng.standard_normal(ha.size))
            for key, val in (("t", t_mjd), ("ha", ha), ("u", u), ("v", v), ("w", w),
                             ("b", np.full(ha.size, b_idx)), ("vis", V),
                             ("wt", np.full(ha.size, 1.0 / sig ** 2))):
                rows[key].append(val)
        cat = lambda k: np.concatenate(rows[k]) if rows[k] else np.zeros(0)
        return VisData(time_mjd=cat("t"), ha_deg=cat("ha"), u=cat("u"), v=cat("v"), w=cat("w"),
                       baseline=cat("b").astype(int), vis=cat("vis"), weight=cat("wt"),
                       freq_hz=[f], baseline_names=obs.array.baseline_names(),
                       header={"OBJECT": obs.name or sky.name, "RA": obs.ra_deg, "DEC": obs.dec_deg,
                               "TELESCOP": obs.array.name, "INSTRUME": rec.name, "BANDWID": bw,
                               "INTTIME": obs.t_int_s, "SIMULATD": True, "NOISE": self.add_noise,
                               "ORIGIN": "interferosim.observation.VisibilitySimulator"})


# ======================================================================
@dataclass
class TransitSimulator:
    """Drift scan of fixed antennas (notebook B, aas6sim ``solar_transit``).

    All antennas point at the meridian at the declination ``dec_deg``
    (fixed in the earth frame). The source -- the field centre of the sky
    model -- moves through the beams; without fringe stopping the
    visibility is sum_c sqrt(P_i P_j) V_c(u, v) exp(-2 pi i w), its phase
    runs with the natural fringe rate. Noise per dump as in
    :class:`VisibilitySimulator`, with the sky power of the moment.
    """

    array: Array
    dec_deg: float
    ha_range_h: Tuple[float, float] = (-2.0, 2.0)
    t_dump_s: float = 5.0
    pointing_ha_deg: float = 0.0
    add_noise: bool = True
    seed: int = 7
    date: str = "2026-11-01"

    def run(self, sky: SkyModel, freq_hz: float = 1.420405751768e9,
            receiver: Optional[Receiver] = None) -> VisData:
        rec = receiver or self.array.receiver
        if rec is None:
            raise ValueError("no receiver")
        if sky.center.ra_deg is None or sky.center.dec_deg is None:
            raise ValueError("the sky model needs a field centre with RA/Dec "
                             "(for the Sun: observation.sun.sun_radec_deg)")
        obs = Observation(self.array, sky.center.ra_deg, sky.center.dec_deg, freq_hz=freq_hz,
                          ha_range_h=self.ha_range_h, t_int_s=self.t_dump_s,
                          min_elev_deg=-90.0, date=self.date, name=sky.center.label or sky.name)
        st = self.array.stations
        ha = obs.hour_angles()
        g = obs.gha_deg(ha)
        # fixed pointing (earth frame): direction of (dec_deg, pointing_ha) at the reference
        p = direction(self.dec_deg, float(obs.gha_deg(np.array([self.pointing_ha_deg]))[0]))
        s0, e_e, e_n = sky_basis(sky.center.dec_deg, g)
        q = obs.uvw(ha)
        sefd = [rec.sefd_jy(s.antenna, freq_hz) if s.antenna else 0.0 for s in st]
        amp = np.zeros((len(st), len(sky), ha.size))
        for ci, c in enumerate(sky):
            l, m = math.radians(c.dx_deg), math.radians(c.dy_deg)
            sc = math.sqrt(max(1 - l * l - m * m, 0)) * s0 + l * e_e + m * e_n
            theta = np.degrees(np.arccos(np.clip(sc @ p, -1.0, 1.0)))
            for k, s in enumerate(st):
                amp[k, ci] = _beam_amp(s.antenna, theta, freq_hz)
        s_in = np.einsum("c,kct->kt", np.array([c.flux_jy for c in sky]), amp ** 2)
        rng = np.random.default_rng(self.seed)
        rows = {k: [] for k in ("t", "ha", "u", "v", "w", "b", "vis", "wt")}
        t_mjd = obs.times_mjd()
        bw = rec.bandwidth_hz
        for b_idx, (i, j) in enumerate(self.array.baseline_pairs()):
            u, v, w = q[b_idx, :, 0], q[b_idx, :, 1], q[b_idx, :, 2]
            V = np.zeros(ha.size, complex)
            for ci, c in enumerate(sky):
                V += amp[i, ci] * amp[j, ci] * c.visibility(u, v)
            V = V * np.exp(-2j * np.pi * w)                      # no fringe stopping
            sig = np.sqrt((sefd[i] + s_in[i]) * (sefd[j] + s_in[j]) / (2.0 * bw * self.t_dump_s))
            if self.add_noise:
                V = V + sig * (rng.standard_normal(ha.size) + 1j * rng.standard_normal(ha.size))
            for key, val in (("t", t_mjd), ("ha", ha), ("u", u), ("v", v), ("w", w),
                             ("b", np.full(ha.size, b_idx)), ("vis", V), ("wt", 1.0 / sig ** 2)):
                rows[key].append(val)
        cat = lambda k: np.concatenate(rows[k])
        return VisData(time_mjd=cat("t"), ha_deg=cat("ha"), u=cat("u"), v=cat("v"), w=cat("w"),
                       baseline=cat("b").astype(int), vis=cat("vis"), weight=cat("wt"),
                       freq_hz=[freq_hz], baseline_names=self.array.baseline_names(),
                       header={"OBJECT": sky.center.label or sky.name, "RA": sky.center.ra_deg,
                               "DEC": sky.center.dec_deg, "POINTDEC": self.dec_deg,
                               "TELESCOP": self.array.name, "INSTRUME": rec.name, "BANDWID": bw,
                               "INTTIME": self.t_dump_s, "SIMULATD": True, "FRINGSTP": False,
                               "ORIGIN": "interferosim.observation.TransitSimulator"})
