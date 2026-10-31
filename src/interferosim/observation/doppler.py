# src/interferosim/observation/doppler.py -- Doppler correction to the LSR (interferosim.observation.doppler)
"""Velocity of the observer and the correction from topocentric to LSR
velocities (notebook 09).

A spectral line is catalogued at its velocity relative to the **local
standard of rest** (LSR). The observer moves relative to the LSR with

* the Earth's rotation (up to 0.46 km/s at the equator, 0.31 km/s at
  Wilhelmsburg), period one sidereal day,
* the Earth's orbit around the Sun (about 30 km/s), period one year,
* the motion of the Sun relative to the LSR (standard solar motion
  20 km/s towards RA 18h03m50.29s, Dec +30d00m16.8s, J2000).

With the unit vector s towards the source, the radial velocities (positive
away from the observer) are related by

    v_LSR = v_topo + v_obs . s ,

where v_obs is the velocity of the observer relative to the LSR
(:func:`v_lsr_correction_kms` returns v_obs . s). The topocentric frequency of
a line with the LSR velocity v (radio convention) is
f = f_0 (1 - (v - v_obs . s) / c) (:func:`topocentric_frequency_hz`).
:func:`regrid_to_lsr` interpolates multichannel visibilities recorded with a
fixed LO onto a common LSR velocity axis (Doppler tracking in software).

Accuracy: orbital velocity from the low-precision solar position (formula of
:mod:`interferosim.observation.sun`, referred to the J2000 equinox,
derivative over +-0.5 d); compared with astropy (barycentric correction plus
the same solar motion) within 0.02 km/s for 2026-2027 -- enough to place a
line in a channel of a few 100 m/s.

Version 0.9.0 (2026-11-01): first implementation.
"""

from __future__ import annotations

import math
from typing import Iterable, Tuple

import numpy as np

from ..constants import C_LIGHT
from ..instrument.array import geodetic_to_ecef
from .geometry import gmst_deg

__all__ = ["observer_velocity_kms", "v_lsr_correction_kms", "topocentric_frequency_hz",
           "lsr_velocity_from_topocentric", "regrid_to_lsr", "SOLAR_MOTION"]

C_KMS = C_LIGHT / 1000.0
AU_KM = 149_597_870.7
OMEGA_EARTH = 7.2921150e-5                         # rad/s
# standard solar motion relative to the LSR: (speed km/s, RA deg, Dec deg), J2000
SOLAR_MOTION = (20.0, 270.9595417, 30.0046667)


def _unit(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec) * np.ones_like(ra)],
                    axis=-1)


def _sun_vector_au(mjd) -> np.ndarray:
    """Geocentric position of the Sun in equatorial coordinates [AU]."""
    n = np.asarray(mjd, dtype=float) - 51544.5
    L = np.radians(np.mod(280.460 + 0.9856474 * n, 360.0))
    g = np.radians(np.mod(357.528 + 0.9856003 * n, 360.0))
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    lam = lam - np.radians(1.3969713 * n / 36525.0)          # ecliptic of date -> J2000 (precession)
    eps = np.radians(23.4392911)                              # J2000 obliquity
    R = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    return np.stack([R * np.cos(lam), R * np.cos(eps) * np.sin(lam), R * np.sin(eps) * np.sin(lam)],
                    axis=-1)


def observer_velocity_kms(mjd, lat_deg: float, lon_deg: float, h_m: float = 0.0,
                          components: Iterable[str] = ("rotation", "orbit", "sun")) -> np.ndarray:
    """Velocity of the observer relative to the LSR in equatorial (J2000-ish)
    coordinates [km/s], shape (3,) or (n, 3). ``components`` selects the
    terms "rotation", "orbit" and "sun" (solar motion)."""
    t = np.asarray(mjd, dtype=float)
    comps = set(components)
    v = np.zeros(t.shape + (3,))
    if "orbit" in comps:
        # Earth relative to the Sun: -d/dt of the geocentric solar vector
        v = v - (_sun_vector_au(t + 0.5) - _sun_vector_au(t - 0.5)) * AU_KM / 86400.0
    if "rotation" in comps:
        r = geodetic_to_ecef(lat_deg, lon_deg, h_m) / 1000.0              # km, ECEF
        th = np.radians(gmst_deg((t - 40587.0) * 86400.0))
        x = r[0] * np.cos(th) - r[1] * np.sin(th)                          # equatorial (inertial)
        y = r[0] * np.sin(th) + r[1] * np.cos(th)
        v = v + np.stack([-OMEGA_EARTH * y, OMEGA_EARTH * x, np.zeros_like(x)], axis=-1)
    if "sun" in comps:
        sp, ra, dec = SOLAR_MOTION
        v = v + sp * _unit(ra, dec)
    return v


def v_lsr_correction_kms(ra_deg: float, dec_deg: float, mjd, lat_deg: float, lon_deg: float,
                         h_m: float = 0.0,
                         components: Iterable[str] = ("rotation", "orbit", "sun")) -> np.ndarray:
    """v_obs . s [km/s]: add to a topocentric radial velocity to get the LSR
    velocity of the source at (RA, Dec)."""
    s = _unit(ra_deg, dec_deg)
    return observer_velocity_kms(mjd, lat_deg, lon_deg, h_m, components) @ s


def topocentric_frequency_hz(rest_hz: float, v_lsr_kms, ra_deg: float, dec_deg: float, mjd,
                             lat_deg: float, lon_deg: float, h_m: float = 0.0) -> np.ndarray:
    """Observed (topocentric) frequency of a line with the LSR velocity
    ``v_lsr_kms`` (radio convention)."""
    vc = v_lsr_correction_kms(ra_deg, dec_deg, mjd, lat_deg, lon_deg, h_m)
    return rest_hz * (1.0 - (np.asarray(v_lsr_kms, dtype=float) - vc) / C_KMS)


def lsr_velocity_from_topocentric(freq_hz, rest_hz: float, ra_deg: float, dec_deg: float, mjd,
                                  lat_deg: float, lon_deg: float, h_m: float = 0.0) -> np.ndarray:
    """LSR radio velocity [km/s] of the topocentric frequencies ``freq_hz``."""
    v_topo = C_KMS * (1.0 - np.asarray(freq_hz, dtype=float) / rest_hz)
    return v_topo + v_lsr_correction_kms(ra_deg, dec_deg, mjd, lat_deg, lon_deg, h_m)


def regrid_to_lsr(vd, array, rest_hz: float, v_axis_kms) -> "VisData":
    """Doppler tracking in software: interpolate every row of ``vd``
    (topocentric channels) onto the common LSR velocity axis ``v_axis_kms``
    (radio convention, rest frequency ``rest_hz``), using the correction for
    the time of the row at the reference station of ``array`` and the field
    centre in ``vd.header`` (RA, DEC). Channels outside the observed band get
    weight 0. The new frequency axis is f = f0 (1 - v/c); u, v, w are
    rescaled to its first channel."""
    from ..visdata import VisData                         # local import (no cycle at import time)
    v_new = np.asarray(v_axis_kms, dtype=float)
    f_new = rest_hz * (1.0 - v_new / C_KMS)
    ref = array.reference
    vc = v_lsr_correction_kms(float(vd.header["RA"]), float(vd.header["DEC"]), vd.time_mjd,
                              ref.lat_deg, ref.lon_deg, ref.height_m)
    vis = np.zeros((vd.n_rows, v_new.size), complex)
    wt = np.zeros((vd.n_rows, v_new.size))
    order = np.argsort(vd.freq_hz)
    f_old = vd.freq_hz[order]
    for r in range(vd.n_rows):
        v_lsr_old = C_KMS * (1.0 - f_old / rest_hz) + vc[r]               # LSR velocity of each channel
        idx = np.argsort(v_lsr_old)
        x = v_lsr_old[idx]
        inside = (v_new >= x[0]) & (v_new <= x[-1])
        row = vd.vis[r, order][idx]
        vis[r] = np.interp(v_new, x, row.real) + 1j * np.interp(v_new, x, row.imag)
        wt[r] = np.where(inside, np.interp(v_new, x, vd.weight[r, order][idx]), 0.0)
    s = f_new[0] / vd.freq_hz[0]
    hdr = dict(vd.header, VELREF="LSR", RESTFREQ=rest_hz, REGRID=True)
    return VisData(time_mjd=vd.time_mjd, ha_deg=vd.ha_deg, u=vd.u * s, v=vd.v * s, w=vd.w * s,
                   baseline=vd.baseline, vis=vis, weight=wt, freq_hz=f_new,
                   baseline_names=list(vd.baseline_names), header=hdr)
