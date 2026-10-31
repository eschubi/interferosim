# src/interferosim/observation/spectral.py -- multichannel visibilities with spectral lines (interferosim.observation.spectral)
"""Visibilities over many frequency channels: continuum and spectral lines
(notebook 09).

:class:`SpectralVisibilitySimulator` extends the single-frequency
:class:`~interferosim.observation.VisibilitySimulator` by a channel axis:

* the band of ``bandwidth_hz`` is split into ``n_chan`` channels with the
  frequencies f_k = f_c + (k - n/2) df (layout of the FX correlator,
  ``Accumulator.freq_axis``);
* the band centre f_c is either given or set from a rest frequency and an
  LSR velocity at the time of transit (the LO stays fixed during the
  observation, as with a real B210 recording);
* every channel is evaluated at its own spatial frequency (u f_k / f_0) and
  with the spectra of the components (:attr:`SkyModel.spectra`); with
  ``doppler`` the topocentric channel frequency is converted to the LSR frame
  for every integration, f_LSR = f_topo (1 - v_obs . s / c), i.e.
  v_LSR = v_topo + v_obs . s, so the lines
  drift through the channels over the day and the year;
* noise per channel and integration
  sqrt((SEFD_i + S_i(f)) (SEFD_j + S_j(f)) / (2 df t_int)).

Result: :class:`~interferosim.visdata.VisData` with ``n_chan`` channels,
u, v, w at ``freq_hz[0]``, header keys RESTFREQ, VLSRCORR (mean correction),
CHANWID.

Version 0.9.0 (2026-11-01): first implementation.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from ..constants import C_LIGHT
from ..instrument import Receiver
from ..sky import SkyModel
from ..visdata import VisData
from .doppler import topocentric_frequency_hz, v_lsr_correction_kms
from .obs import Observation, _beam_amp

__all__ = ["SpectralVisibilitySimulator"]

C_KMS = C_LIGHT / 1000.0


class SpectralVisibilitySimulator:
    """Analytic multichannel visibilities of a sky model with spectra."""

    def __init__(self, receiver: Optional[Receiver] = None, *, n_chan: int = 64,
                 bandwidth_hz: Optional[float] = None, center_freq_hz: Optional[float] = None,
                 rest_hz: Optional[float] = None, v_center_kms: float = 0.0,
                 doppler: bool = True, add_noise: bool = True, seed: int = 1):
        if n_chan < 1:
            raise ValueError("n_chan must be >= 1")
        self.receiver = receiver
        self.n_chan = int(n_chan)
        self.bandwidth_hz = bandwidth_hz
        self.center_freq_hz = center_freq_hz
        self.rest_hz = rest_hz
        self.v_center_kms = v_center_kms
        self.doppler = doppler
        self.add_noise = add_noise
        self.seed = seed

    # ------------------------------------------------------------------
    def band_center_hz(self, sky: SkyModel, obs: Observation) -> float:
        """Fixed LO: given centre, or the topocentric frequency of
        ``v_center_kms`` at ``rest_hz`` at the time of transit."""
        if self.center_freq_hz:
            return float(self.center_freq_hz)
        rest = self.rest_hz or obs.freq_hz
        if not self.doppler:
            return float(rest * (1.0 - self.v_center_kms / C_KMS))
        ref = obs.array.reference
        return float(topocentric_frequency_hz(rest, self.v_center_kms, sky.center.ra_deg,
                                              sky.center.dec_deg, obs.transit_mjd(),
                                              ref.lat_deg, ref.lon_deg, ref.height_m))

    def freq_axis(self, sky: SkyModel, obs: Observation, bandwidth_hz: float) -> np.ndarray:
        df = bandwidth_hz / self.n_chan
        return self.band_center_hz(sky, obs) + (np.arange(self.n_chan) - self.n_chan // 2) * df

    # ------------------------------------------------------------------
    def run(self, sky: SkyModel, obs: Observation, only_visible: bool = True) -> VisData:
        rec = self.receiver or obs.array.receiver
        if rec is None:
            raise ValueError("no receiver")
        bw = self.bandwidth_hz or rec.bandwidth_hz
        df = bw / self.n_chan
        freq = self.freq_axis(sky, obs, bw)
        st = obs.array.stations
        ha = obs.hour_angles()
        keep = obs.visible(ha) if only_visible else np.ones(ha.size, bool)
        ha = ha[keep]
        t_mjd = obs.times_mjd()[keep]
        q = obs.uvw(ha) * (freq[0] / obs.freq_hz)                          # at freq_hz[0]
        ref = obs.array.reference
        if self.doppler and ha.size:
            vc = v_lsr_correction_kms(sky.center.ra_deg, sky.center.dec_deg, t_mjd,
                                      ref.lat_deg, ref.lon_deg, ref.height_m)
        else:
            vc = np.zeros(ha.size)
        f_lsr = freq[None, :] * (1.0 - vc[:, None] / C_KMS)                # v_LSR = v_topo + vc
        # shape of every component per integration and channel
        shapes = []
        for c in sky:
            sp = sky.spectra.get(c.name)
            shapes.append(sp.shape(f_lsr) if sp is not None else np.ones_like(f_lsr))
        theta = [math.hypot(c.dx_deg, c.dy_deg) for c in sky]
        amp = np.array([[float(_beam_amp(s.antenna, th, obs.freq_hz)) for th in theta] for s in st])
        s_in = np.zeros((len(st), ha.size, self.n_chan))
        for k in range(len(st)):
            for ci, c in enumerate(sky):
                s_in[k] += c.flux_jy * amp[k, ci] ** 2 * shapes[ci]
        sefd = [rec.sefd_jy(s.antenna, obs.freq_hz) if s.antenna else 0.0 for s in st]
        rng = np.random.default_rng(self.seed)
        rows = {k: [] for k in ("t", "ha", "u", "v", "w", "b", "vis", "wt")}
        for b_idx, (i, j) in enumerate(obs.array.baseline_pairs()):
            u, v, w = q[b_idx, :, 0], q[b_idx, :, 1], q[b_idx, :, 2]
            V = np.zeros((ha.size, self.n_chan), complex)
            for kc, fk in enumerate(freq):
                s = fk / freq[0]
                for ci, c in enumerate(sky):
                    if c.flux_jy == 0:
                        continue
                    V[:, kc] += amp[i, ci] * amp[j, ci] * c.visibility(u * s, v * s) * shapes[ci][:, kc]
            sig = np.sqrt((sefd[i] + s_in[i]) * (sefd[j] + s_in[j]) / (2.0 * df * obs.t_int_s))
            if self.add_noise:
                V = V + sig * (rng.standard_normal(V.shape) + 1j * rng.standard_normal(V.shape))
            for key, val in (("t", t_mjd), ("ha", ha), ("u", u), ("v", v), ("w", w),
                             ("b", np.full(ha.size, b_idx)), ("vis", V), ("wt", 1.0 / sig ** 2)):
                rows[key].append(val)
        cat = lambda k: np.concatenate(rows[k]) if rows[k] else np.zeros(0)
        vis = np.concatenate(rows["vis"]) if rows["vis"] else np.zeros((0, self.n_chan), complex)
        wt = np.concatenate(rows["wt"]) if rows["wt"] else np.zeros((0, self.n_chan))
        return VisData(time_mjd=cat("t"), ha_deg=cat("ha"), u=cat("u"), v=cat("v"), w=cat("w"),
                       baseline=cat("b").astype(int), vis=vis, weight=wt, freq_hz=freq,
                       baseline_names=obs.array.baseline_names(),
                       header={"OBJECT": obs.name or sky.name, "RA": sky.center.ra_deg,
                               "DEC": sky.center.dec_deg, "TELESCOP": obs.array.name,
                               "INSTRUME": rec.name, "BANDWID": bw, "CHANWID": df,
                               "RESTFREQ": self.rest_hz or obs.freq_hz,
                               "VLSRCORR": float(np.mean(vc)) if vc.size else 0.0,
                               "DOPPLER": self.doppler, "INTTIME": obs.t_int_s, "SIMULATD": True,
                               "NOISE": self.add_noise,
                               "ORIGIN": "interferosim.observation.SpectralVisibilitySimulator"})
