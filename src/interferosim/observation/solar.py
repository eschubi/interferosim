# src/interferosim/observation/solar.py -- the Sun as a target: flux, transit, disc fit (interferosim.observation.solar)
"""The Sun as the first target of the project (notebook 08, Appendix B of
the teaching material).

* :func:`sun_flux_jy` -- flux density of a disc of brightness temperature
  T_b and diameter theta: S = 2 k T_b Omega / lambda^2.
* :func:`solar_transit_mjd` -- time of the upper transit of the Sun (local
  apparent noon) at a longitude on a date.
* :class:`SolarTransitSimulator` -- drift scan of fixed antennas through
  which the Sun moves. Unlike :class:`~interferosim.observation.TransitSimulator`
  (fixed RA/Dec, sidereal rate) the position of the Sun is recomputed for
  every dump: its hour angle advances with the *solar* day, 0.27 % slower
  than a fixed star -- over two hours of a 10 m baseline about 0.07 fringe.
* :func:`disc_amplitude` and :func:`fit_disc_diameter` -- visibility
  amplitude of a uniform disc, ``|2 J1(x)/x|`` with x = pi theta q, and its
  least-squares fit to measured amplitudes over the baseline length
  (diameter measurement of notebook B).

Version 0.8.1 (2026-10-31): first implementation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import j1

from ..constants import C_LIGHT, JANSKY, K_BOLTZMANN
from ..instrument import Array, Receiver
from ..sky import SkyModel
from ..visdata import VisData
from .geometry import direction, gmst_deg, sky_basis, uvw
from .obs import _beam_amp
from .sun import mjd_from_iso, sun_radec_deg

__all__ = ["sun_flux_jy", "solar_transit_mjd", "SolarTransitSimulator", "disc_amplitude",
           "fit_disc_diameter"]


def sun_flux_jy(t_b_k: float, diameter_deg: float, freq_hz: float = 1.420405751768e9) -> float:
    """Flux density [Jy] of a uniform disc, S = 2 k T_b Omega / lambda^2
    (Rayleigh-Jeans), Omega = pi (theta/2)^2."""
    lam = C_LIGHT / freq_hz
    omega = math.pi * (0.5 * math.radians(diameter_deg)) ** 2
    return 2.0 * K_BOLTZMANN * t_b_k * omega / lam ** 2 / JANSKY


def _unix(mjd):
    return (np.asarray(mjd, dtype=float) - 40587.0) * 86400.0


def _lha_sun_deg(mjd, lon_deg: float) -> np.ndarray:
    ra, _ = sun_radec_deg(mjd)
    return np.mod(gmst_deg(_unix(mjd)) + lon_deg - ra + 180.0, 360.0) - 180.0


def solar_transit_mjd(date: str, lon_deg: float) -> float:
    """MJD (UTC) of the upper transit of the Sun at ``lon_deg`` (east
    positive) on the UTC date ``date`` (local apparent noon)."""
    mjd = mjd_from_iso(date) + 0.5 - lon_deg / 360.0
    for _ in range(4):
        mjd -= float(_lha_sun_deg(mjd, lon_deg)) / 360.0 * (1.0 - 1.0 / 365.2422)
    return mjd


# ======================================================================
@dataclass
class SolarTransitSimulator:
    """Fixed antennas, the Sun drifts through the beams (notebook 08).

    All antennas point at the hour angle ``pointing_ha_deg`` (default: the
    meridian) and the declination ``pointing_dec_deg`` (default: that of the
    Sun at transit), fixed in the earth frame. ``ha_range_h`` is given in
    hours of solar time around the transit. For every dump the RA/Dec of the
    Sun are recomputed; the sky model is placed at that position (offsets of
    the components relative to the centre of the disc). Visibility without
    fringe stopping: sum_c sqrt(P_i P_j) V_c(u, v) exp(-2 pi i w); noise per
    dump with the sky power of the moment.
    """

    array: Array
    date: str = "2026-11-01"
    ha_range_h: Tuple[float, float] = (-2.0, 2.0)
    t_dump_s: float = 5.0
    pointing_ha_deg: float = 0.0
    pointing_dec_deg: Optional[float] = None
    add_noise: bool = True
    seed: int = 7

    def transit_mjd(self) -> float:
        return solar_transit_mjd(self.date, self.array.reference.lon_deg)

    def times_mjd(self) -> np.ndarray:
        h0, h1 = self.ha_range_h
        n = max(int(round((h1 - h0) * 3600.0 / self.t_dump_s)), 1)
        return self.transit_mjd() + (h0 * 3600.0 + (np.arange(n) + 0.5) * self.t_dump_s) / 86400.0

    def sun_track(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """(t_mjd, lha_deg, ra_deg, dec_deg) of the Sun for every dump."""
        t = self.times_mjd()
        ra, dec = sun_radec_deg(t)
        lon = self.array.reference.lon_deg
        lha = np.mod(gmst_deg(_unix(t)) + lon - ra + 180.0, 360.0) - 180.0
        return t, lha, ra, dec

    def run(self, sky: SkyModel, freq_hz: float = 1.420405751768e9,
            receiver: Optional[Receiver] = None) -> VisData:
        rec = receiver or self.array.receiver
        if rec is None:
            raise ValueError("no receiver")
        st = self.array.stations
        lon = self.array.reference.lon_deg
        lam = C_LIGHT / freq_hz
        t, lha, ra, dec = self.sun_track()
        g = np.mod(lha - lon, 360.0)                                   # Greenwich hour angle
        dec0 = float(np.mean(dec))                                     # changes < 0.1 deg in hours
        pdec = dec0 if self.pointing_dec_deg is None else self.pointing_dec_deg
        p = direction(pdec, (self.pointing_ha_deg - lon) % 360.0)      # fixed pointing (earth frame)
        s0, e_e, e_n = sky_basis(dec0, g)
        amp = np.zeros((len(st), len(sky), t.size))
        for ci, c in enumerate(sky):
            l, m = math.radians(c.dx_deg), math.radians(c.dy_deg)
            sc = math.sqrt(max(1 - l * l - m * m, 0.0)) * s0 + l * e_e + m * e_n
            theta = np.degrees(np.arccos(np.clip(sc @ p, -1.0, 1.0)))
            for k, s in enumerate(st):
                amp[k, ci] = _beam_amp(s.antenna, theta, freq_hz)
        s_in = np.einsum("c,kct->kt", np.array([c.flux_jy for c in sky]), amp ** 2)
        sefd = [rec.sefd_jy(s.antenna, freq_hz) if s.antenna else 0.0 for s in st]
        rng = np.random.default_rng(self.seed)
        bw = rec.bandwidth_hz
        rows = {k: [] for k in ("t", "ha", "u", "v", "w", "b", "vis", "wt")}
        for b_idx, (i, j) in enumerate(self.array.baseline_pairs()):
            q = uvw(self.array.baselines()[b_idx], dec0, g, lam)
            u, v, w = q[:, 0], q[:, 1], q[:, 2]
            V = np.zeros(t.size, complex)
            for ci, c in enumerate(sky):
                V += amp[i, ci] * amp[j, ci] * c.visibility(u, v)
            V = V * np.exp(-2j * np.pi * w)                              # no fringe stopping
            sig = np.sqrt((sefd[i] + s_in[i]) * (sefd[j] + s_in[j]) / (2.0 * bw * self.t_dump_s))
            if self.add_noise:
                V = V + sig * (rng.standard_normal(t.size) + 1j * rng.standard_normal(t.size))
            for key, val in (("t", t), ("ha", lha), ("u", u), ("v", v), ("w", w),
                             ("b", np.full(t.size, b_idx)), ("vis", V), ("wt", 1.0 / sig ** 2)):
                rows[key].append(val)
        cat = lambda k: np.concatenate(rows[k])
        return VisData(time_mjd=cat("t"), ha_deg=cat("ha"), u=cat("u"), v=cat("v"), w=cat("w"),
                       baseline=cat("b").astype(int), vis=cat("vis"), weight=cat("wt"),
                       freq_hz=[freq_hz], baseline_names=self.array.baseline_names(),
                       header={"OBJECT": "Sun", "RA": float(ra[ra.size // 2]), "DEC": dec0,
                               "POINTDEC": pdec, "TELESCOP": self.array.name, "INSTRUME": rec.name,
                               "BANDWID": bw, "INTTIME": self.t_dump_s, "SIMULATD": True,
                               "FRINGSTP": False, "SOLARDAY": True,
                               "ORIGIN": "interferosim.observation.SolarTransitSimulator"})


# ======================================================================
def disc_amplitude(q_lambda, s0: float, diameter_deg: float) -> np.ndarray:
    """``|V|`` of a uniform disc: ``s0 |2 J1(x)/x|``, x = pi theta q (q in wavelengths)."""
    x = math.pi * math.radians(diameter_deg) * np.asarray(q_lambda, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return s0 * np.abs(np.where(x < 1e-9, 1.0, 2.0 * j1(x) / np.where(x == 0, 1.0, x)))


def fit_disc_diameter(q_lambda: Sequence[float], amplitude: Sequence[float],
                      sigma: Optional[Sequence[float]] = None,
                      p0: Tuple[float, float] = (1.0, 0.5)) -> Tuple[float, float, float, float]:
    """Least-squares fit of :func:`disc_amplitude` to measured amplitudes.
    Returns (s0, diameter_deg, error of s0, error of the diameter)."""
    q = np.asarray(q_lambda, dtype=float)
    a = np.asarray(amplitude, dtype=float)
    p, cov = curve_fit(lambda x, s, d: disc_amplitude(x, s, d), q, a, p0=p0,
                       sigma=None if sigma is None else np.asarray(sigma, dtype=float),
                       absolute_sigma=sigma is not None)
    e = np.sqrt(np.diag(cov))
    return float(p[0]), float(p[1]), float(e[0]), float(e[1])
