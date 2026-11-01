# src/interferosim/observation/geometry.py -- source direction, geometric delay, uvw (interferosim.observation.geometry)
"""Time-dependent geometry of an interferometer baseline.

The source direction is computed in the earth-fixed frame (ECEF) from its
right ascension and declination and the Greenwich hour angle
GHA = GMST - RA. All baselines are ECEF vectors (``Array.baselines()``).

* :func:`gmst_deg` -- Greenwich mean sidereal time from UTC (unix time).
* :func:`sky_basis` -- unit vectors towards the phase centre (s0) and east
  (e_E, increasing RA) and north (e_N) on the sky, in ECEF.
* :func:`direction` -- unit vector of a source at the offset (l, m) from the
  phase centre (small-field: l = dx, m = dy in radians).
* :func:`geometric_delay` -- tau = -(B . s) / c.
* :func:`uvw` -- baseline in wavelengths, projected on (e_E, e_N, s0).
* :func:`elevation_deg` -- elevation of a direction at a site.

Added in 0.9.4 (theory notebook E):

* :func:`lst_deg` -- local sidereal time, GMST + longitude.
* :func:`azimuth_elevation_deg` -- azimuth (from north through east) and
  elevation of a direction at a site.
* :func:`parallactic_angle_deg` -- angle between the direction to the
  zenith and to the north celestial pole at the source,
  q = atan2(sin H, tan(phi) cos(delta) - sin(delta) cos H).
* :func:`precess_deg` -- precession of J2000 coordinates to the mean
  equator and equinox of a date (IAU 1976 angles zeta, z, theta). The
  simulators use the catalogue coordinates directly with the GMST of the
  date (self-consistent: the correlator uses the same model); for real data
  the coordinates must be referred to the date (notebook E shows the size).

Conventions (as AAS6 ``modules/interfero/geometry.py``)
-------------------------------------------------------
* B = r_j - r_i; station i is the reference.
* tau > 0: the wavefront reaches station j LATER than station i.
* For a baseband signal with the LO at f_c the delay appears as the phase
  factor exp(-2 pi i f_c tau) = exp(+2 pi i (B . s) / lambda); the cross
  product x_i x_j* of a source at (l, m) is therefore
  S exp(-2 pi i (w + u l + v m)), and after removing the phase-centre term
  (fringe stopping) S exp(-2 pi i (u l + v m)) -- the convention of
  ``SkyModel.visibility``.

AAS6 uses astropy (ITRS) for the source direction; here the classical
GMST formula is used (IAU 1982, accuracy ~0.1 s over decades), without
precession/nutation -- sufficient for the simulation, where the correlator
uses the same model.

Version 0.4.0 (2026-10-31): first implementation.
Version 0.9.4 (2026-11-01): lst_deg, azimuth_elevation_deg,
parallactic_angle_deg, precess_deg.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

from ..constants import C_LIGHT, OMEGA_SIDEREAL

__all__ = ["gmst_deg", "gha_deg", "sky_basis", "direction", "geometric_delay",
           "delay_rate", "uvw", "elevation_deg", "SIDEREAL_DEG_PER_S",
           "lst_deg", "azimuth_elevation_deg", "parallactic_angle_deg", "precess_deg"]

SIDEREAL_DEG_PER_S = math.degrees(OMEGA_SIDEREAL)     # 0.004178 deg/s


def gmst_deg(unix_time_s) -> np.ndarray:
    """Greenwich mean sidereal time in degrees for UTC given as unix time."""
    jd = np.asarray(unix_time_s, dtype=float) / 86400.0 + 2440587.5
    d = jd - 2451545.0
    return np.mod(280.46061837 + 360.98564736629 * d, 360.0)


def gha_deg(ra_deg: float, unix_time_s=None, *, lha_deg=None, lon_deg: float = 0.0):
    """Greenwich hour angle of a source.

    Either from UTC (``unix_time_s``) or from the local hour angle
    ``lha_deg`` at a site of longitude ``lon_deg`` (east positive):
    GHA = LHA - lon.
    """
    if lha_deg is not None:
        return np.mod(np.asarray(lha_deg, dtype=float) - lon_deg, 360.0)
    return np.mod(gmst_deg(unix_time_s) - ra_deg, 360.0)


def sky_basis(dec_deg: float, gha) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit vectors (s0, e_E, e_N) in ECEF for declination ``dec_deg`` and
    Greenwich hour angle ``gha`` (degrees, scalar or array of length n).
    Returns arrays of shape (3,) or (n, 3)."""
    h = np.radians(np.asarray(gha, dtype=float))
    d = math.radians(dec_deg)
    ch, sh = np.cos(h), np.sin(h)
    s0 = np.stack([math.cos(d) * ch, -math.cos(d) * sh, np.full_like(h, math.sin(d))], axis=-1)
    e_e = np.stack([sh, ch, np.zeros_like(h)], axis=-1)
    e_n = np.stack([-math.sin(d) * ch, math.sin(d) * sh, np.full_like(h, math.cos(d))], axis=-1)
    return s0, e_e, e_n


def direction(dec_deg: float, gha, l_rad: float = 0.0, m_rad: float = 0.0) -> np.ndarray:
    """Unit vector towards a source at the offset (l, m) [rad] (east,
    north) from the phase centre (dec_deg, gha)."""
    s0, e_e, e_n = sky_basis(dec_deg, gha)
    n_ = math.sqrt(max(1.0 - l_rad * l_rad - m_rad * m_rad, 0.0))
    return n_ * s0 + l_rad * e_e + m_rad * e_n


def geometric_delay(baseline_ecef_m, s) -> np.ndarray:
    """Geometric delay tau = -(B . s) / c in seconds (s: unit vector(s))."""
    return -np.asarray(s) @ np.asarray(baseline_ecef_m, dtype=float) / C_LIGHT


def delay_rate(baseline_ecef_m, dec_deg: float, gha, l_rad: float = 0.0,
               m_rad: float = 0.0) -> np.ndarray:
    """Time derivative d tau / dt (s/s) from the Earth's rotation."""
    eps = 1e-3                                         # degrees of GHA
    g = np.asarray(gha, dtype=float)
    t1 = geometric_delay(baseline_ecef_m, direction(dec_deg, g + eps, l_rad, m_rad))
    t0 = geometric_delay(baseline_ecef_m, direction(dec_deg, g - eps, l_rad, m_rad))
    return (t1 - t0) / (2 * eps / SIDEREAL_DEG_PER_S)


def uvw(baseline_ecef_m, dec_deg: float, gha, wavelength_m: float) -> np.ndarray:
    """(u, v, w) of a baseline in wavelengths, shape (3,) or (n, 3)."""
    s0, e_e, e_n = sky_basis(dec_deg, gha)
    b = np.asarray(baseline_ecef_m, dtype=float)
    return np.stack([e_e @ b, e_n @ b, s0 @ b], axis=-1) / wavelength_m


def elevation_deg(s, lat_deg: float, lon_deg: float) -> np.ndarray:
    """Elevation of the direction ``s`` (ECEF unit vector(s)) at a site."""
    p, l = math.radians(lat_deg), math.radians(lon_deg)
    up = np.array([math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p)])
    return np.degrees(np.arcsin(np.clip(np.asarray(s) @ up, -1.0, 1.0)))


# -- added in 0.9.4 (theory notebook E) -----------------------------------------
def lst_deg(unix_time_s, lon_deg: float) -> np.ndarray:
    """Local (mean) sidereal time in degrees: GMST + east longitude."""
    return np.mod(gmst_deg(unix_time_s) + lon_deg, 360.0)


def azimuth_elevation_deg(s, lat_deg: float, lon_deg: float) -> Tuple[np.ndarray, np.ndarray]:
    """(azimuth, elevation) in degrees of the ECEF unit vector(s) ``s`` at a
    site; azimuth from north through east, 0..360."""
    p, l = math.radians(lat_deg), math.radians(lon_deg)
    east = np.array([-math.sin(l), math.cos(l), 0.0])
    north = np.array([-math.sin(p) * math.cos(l), -math.sin(p) * math.sin(l), math.cos(p)])
    up = np.array([math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p)])
    s = np.asarray(s, dtype=float)
    e, n, u = s @ east, s @ north, s @ up
    return np.mod(np.degrees(np.arctan2(e, n)), 360.0), np.degrees(np.arcsin(np.clip(u, -1.0, 1.0)))


def parallactic_angle_deg(ha_deg, dec_deg: float, lat_deg: float) -> np.ndarray:
    """Parallactic angle q in degrees (positive west of the meridian):
    q = atan2(sin H, tan(phi) cos(delta) - sin(delta) cos H)."""
    h = np.radians(np.asarray(ha_deg, dtype=float))
    d, p = math.radians(dec_deg), math.radians(lat_deg)
    return np.degrees(np.arctan2(np.sin(h), math.tan(p) * math.cos(d) - math.sin(d) * np.cos(h)))


def precess_deg(ra_deg, dec_deg, mjd) -> Tuple[np.ndarray, np.ndarray]:
    """Precess J2000 coordinates to the mean equator and equinox of ``mjd``
    (IAU 1976: zeta, z, theta as cubic polynomials in Julian centuries)."""
    t = (np.asarray(mjd, dtype=float) - 51544.5) / 36525.0
    asec = math.pi / (180.0 * 3600.0)
    zeta = (2306.2181 * t + 0.30188 * t ** 2 + 0.017998 * t ** 3) * asec
    z = (2306.2181 * t + 1.09468 * t ** 2 + 0.018203 * t ** 3) * asec
    theta = (2004.3109 * t - 0.42665 * t ** 2 - 0.041833 * t ** 3) * asec
    a, d = np.radians(np.asarray(ra_deg, dtype=float)), np.radians(np.asarray(dec_deg, dtype=float))
    A = np.cos(d) * np.sin(a + zeta)
    B = np.cos(theta) * np.cos(d) * np.cos(a + zeta) - np.sin(theta) * np.sin(d)
    C = np.sin(theta) * np.cos(d) * np.cos(a + zeta) + np.cos(theta) * np.sin(d)
    return np.mod(np.degrees(np.arctan2(A, B) + z), 360.0), np.degrees(np.arcsin(np.clip(C, -1.0, 1.0)))
