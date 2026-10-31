# src/interferosim/observation/sun.py -- position of the Sun (interferosim.observation.sun)
"""Apparent position of the Sun, low precision.

Formula of the Astronomical Almanac ("Low precision formulas for the Sun",
section C), accuracy about 0.01 deg for 1950-2050 -- sufficient to place
the solar disc (34') in the beam of a horn or a small dish.

* :func:`mjd_from_iso` -- MJD (UTC) from an ISO date/time string.
* :func:`sun_radec_deg` -- right ascension and declination of the Sun.

Version 0.5.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Tuple

import numpy as np

__all__ = ["mjd_from_iso", "sun_radec_deg"]


def mjd_from_iso(text: str) -> float:
    """MJD of an ISO date ("2026-11-01") or date and time
    ("2026-11-01T11:30:00"), interpreted as UTC."""
    dt = datetime.fromisoformat(text)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp() / 86400.0 + 40587.0


def sun_radec_deg(mjd) -> Tuple[np.ndarray, np.ndarray]:
    """(RA, Dec) of the Sun in degrees (J2000-ish, apparent, low precision)."""
    n = np.asarray(mjd, dtype=float) - 51544.5            # days from J2000.0
    L = np.radians(np.mod(280.460 + 0.9856474 * n, 360.0))
    g = np.radians(np.mod(357.528 + 0.9856003 * n, 360.0))
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 4.0e-7 * n)
    ra = np.degrees(np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))) % 360.0
    dec = np.degrees(np.arcsin(np.sin(eps) * np.sin(lam)))
    return ra, dec
