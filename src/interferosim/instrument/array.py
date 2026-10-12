# src/interferosim/instrument/array.py -- site, stations, array (interferosim.instrument.array)
"""Site, stations and array of an interferometer.

* :class:`Site` -- geodetic reference position (WGS84), e.g. Wilhelmsburg.
* :class:`Station` -- one antenna at a geodetic position.
* :class:`Array` -- site, stations and receiver; provides the baselines in
  the earth-fixed frame (ECEF) and in the local frame (east, north, up).

Origin
------
WGS84 conversions from AAS6 ``modules/interfero/geometry.py`` 1.0.0
(``geodetic_to_ecef``, ``ecef_to_geodetic``, ``baseline_vector``),
here without astropy. The time-dependent quantities (source direction,
geometric delay, uvw) follow in the package ``observation``.

Configuration (``configs/sites/*.toml``)
---------------------------------------
A station is given either geodetically (``latitude_deg``,
``longitude_deg``, ``height_m``) or as offset from the site in metres
(``east_m``, ``north_m``, ``up_m``). ``site``, ``antenna`` and
``receiver`` name files in ``configs/sites``, ``configs/antennas`` and
``configs/receivers`` (file name without ``.toml``); a station may name
its own ``antenna``.

Conventions
-----------
* Baseline i-j: B = r_j - r_i (the station listed first is the reference,
  as in AAS6: B = r_B - r_A).
* ENU axes at the site: east, north, up.

Version 0.3.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .antennas import Antenna
from .receiver import Receiver

__all__ = ["Site", "Station", "Array", "geodetic_to_ecef", "ecef_to_geodetic",
           "ecef_to_enu_matrix"]

# WGS84 (AAS6 geometry.py)
_A = 6_378_137.0                 # semi-major axis (m)
_F = 1.0 / 298.257223563         # flattening
_E2 = _F * (2.0 - _F)            # first eccentricity squared


def geodetic_to_ecef(lat_deg: float, lon_deg: float, h_m: float) -> np.ndarray:
    """WGS84 (lat, lon in degrees, height in m) -> ECEF [x, y, z] in m."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    s = math.sin(lat)
    n = _A / math.sqrt(1.0 - _E2 * s * s)
    x = (n + h_m) * math.cos(lat) * math.cos(lon)
    y = (n + h_m) * math.cos(lat) * math.sin(lon)
    z = (n * (1.0 - _E2) + h_m) * s
    return np.array([x, y, z], dtype=float)


def ecef_to_geodetic(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """ECEF [m] -> (lat_deg, lon_deg, h_m), Bowring iteration (AAS6)."""
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - _E2))
    for _ in range(6):
        s = math.sin(lat)
        n = _A / math.sqrt(1.0 - _E2 * s * s)
        h = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1.0 - _E2 * n / (n + h)))
    s = math.sin(lat)
    n = _A / math.sqrt(1.0 - _E2 * s * s)
    h = p / math.cos(lat) - n
    return math.degrees(lat), math.degrees(lon), h


def ecef_to_enu_matrix(lat_deg: float, lon_deg: float) -> np.ndarray:
    """Rotation matrix R with enu = R @ ecef_vector (rows: east, north, up)."""
    p, l = math.radians(lat_deg), math.radians(lon_deg)
    return np.array([
        [-math.sin(l), math.cos(l), 0.0],
        [-math.sin(p) * math.cos(l), -math.sin(p) * math.sin(l), math.cos(p)],
        [math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p)],
    ])


# ======================================================================
@dataclass
class Site:
    """Geodetic reference position of an instrument (WGS84)."""

    name: str
    lat_deg: float
    lon_deg: float
    alt_m: float = 0.0

    def ecef_m(self) -> np.ndarray:
        return geodetic_to_ecef(self.lat_deg, self.lon_deg, self.alt_m)

    def enu_to_ecef(self, east_m: float, north_m: float, up_m: float = 0.0) -> np.ndarray:
        """ECEF position of a point given as offset (east, north, up) from the site."""
        R = ecef_to_enu_matrix(self.lat_deg, self.lon_deg)
        return self.ecef_m() + R.T @ np.array([east_m, north_m, up_m], dtype=float)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Site":
        return cls(name=str(d.get("name", "site")),
                   lat_deg=float(d["latitude_deg"]), lon_deg=float(d["longitude_deg"]),
                   alt_m=float(d.get("altitude_m", d.get("height_m", 0.0))))

    @classmethod
    def from_toml(cls, path: Union[str, Path]) -> "Site":
        """Read ``configs/sites/<name>.toml`` (keys latitude_deg,
        longitude_deg, altitude_m)."""
        with open(path, "rb") as fh:
            return cls.from_dict(tomllib.load(fh))


# ======================================================================
@dataclass
class Station:
    """One antenna at a geodetic position."""

    name: str
    lat_deg: float
    lon_deg: float
    height_m: float = 0.0
    antenna: Optional[Antenna] = None

    def ecef_m(self) -> np.ndarray:
        return geodetic_to_ecef(self.lat_deg, self.lon_deg, self.height_m)

    def enu_m(self, site: Site) -> np.ndarray:
        """Offset (east, north, up) from ``site`` in metres."""
        R = ecef_to_enu_matrix(site.lat_deg, site.lon_deg)
        return R @ (self.ecef_m() - site.ecef_m())

    @classmethod
    def from_enu(cls, name: str, site: Site, east_m: float, north_m: float,
                 up_m: float = 0.0, antenna: Optional[Antenna] = None) -> "Station":
        """Station at an offset (east, north, up) from ``site``."""
        lat, lon, h = ecef_to_geodetic(*site.enu_to_ecef(east_m, north_m, up_m))
        return cls(name, lat, lon, h, antenna)


# ======================================================================
@dataclass
class Array:
    """Site, stations and receiver of an interferometer."""

    site: Site
    stations: List[Station] = field(default_factory=list)
    receiver: Optional[Receiver] = None
    name: str = ""

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.stations)

    def __getitem__(self, name: str) -> Station:
        for s in self.stations:
            if s.name == name:
                return s
        raise KeyError(name)

    @property
    def reference(self) -> Station:
        """Reference station (the first one)."""
        return self.stations[0]

    def baseline_pairs(self) -> List[Tuple[int, int]]:
        """Index pairs (i, j), i < j, of all baselines."""
        return list(combinations(range(len(self.stations)), 2))

    def baseline_names(self) -> List[str]:
        return [f"{self.stations[i].name}-{self.stations[j].name}"
                for i, j in self.baseline_pairs()]

    def baselines(self) -> np.ndarray:
        """Baseline vectors B = r_j - r_i in ECEF, shape (n_baselines, 3), m."""
        r = np.array([s.ecef_m() for s in self.stations])
        return np.array([r[j] - r[i] for i, j in self.baseline_pairs()]).reshape(-1, 3)

    def baselines_enu(self) -> np.ndarray:
        """Baseline vectors in the local frame (east, north, up) of the site."""
        R = ecef_to_enu_matrix(self.site.lat_deg, self.site.lon_deg)
        return self.baselines() @ R.T

    def baseline_lengths_m(self) -> np.ndarray:
        return np.linalg.norm(self.baselines(), axis=1)

    def station_enu(self) -> np.ndarray:
        """Station positions (east, north, up) relative to the site, shape (n, 3)."""
        return np.array([s.enu_m(self.site) for s in self.stations]).reshape(-1, 3)

    # ------------------------------------------------------------------
    @classmethod
    def from_toml(cls, path: Union[str, Path], antenna: Optional[Antenna] = None,
                  receiver: Optional[Receiver] = None) -> "Array":
        """Read an array from ``configs/sites/<name>.toml``.

        ``antenna`` / ``receiver`` override the entries of the file. Named
        files are looked up next to ``path`` (sites) and in the sibling
        directories ``antennas`` and ``receivers``.
        """
        path = Path(path)
        with open(path, "rb") as fh:
            d = tomllib.load(fh)
        configs = path.parent.parent

        def _antenna(name: Optional[str]) -> Optional[Antenna]:
            if not name:
                return None
            return Antenna.from_toml(configs / "antennas" / f"{name}.toml")

        if receiver is None and d.get("receiver"):
            receiver = Receiver.from_toml(configs / "receivers" / f"{d['receiver']}.toml")
        default_ant = antenna or _antenna(d.get("antenna"))

        site: Optional[Site] = None
        if d.get("site"):
            site = Site.from_toml(path.parent / f"{d['site']}.toml")

        st_def = d.get("stations", {})
        if not st_def:
            raise ValueError(f"{path}: no [stations.*] entries")
        ref = d.get("reference")
        names = list(st_def)
        if ref in names:                              # reference first
            names.remove(ref)
            names.insert(0, ref)
        stations: List[Station] = []
        for n in names:
            s = st_def[n]
            ant = antenna or _antenna(s.get("antenna")) or default_ant
            if "latitude_deg" in s:
                stations.append(Station(n, float(s["latitude_deg"]), float(s["longitude_deg"]),
                                        float(s.get("height_m", 0.0)), ant))
            else:
                if site is None:
                    raise ValueError(f"{path}: station {n} given in east/north/up "
                                     f"needs a 'site'")
                stations.append(Station.from_enu(n, site, float(s.get("east_m", 0.0)),
                                                 float(s.get("north_m", 0.0)),
                                                 float(s.get("up_m", 0.0)), ant))
        if site is None:                              # site = reference station
            r = stations[0]
            site = Site(r.name, r.lat_deg, r.lon_deg, r.height_m)
        return cls(site=site, stations=stations, receiver=receiver,
                   name=str(d.get("name", path.stem)))
