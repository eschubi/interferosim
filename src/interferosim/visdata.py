# src/interferosim/visdata.py -- visibility data class (interferosim.visdata)
"""Visibility data shared by observation, signal, imaging and io.

:class:`VisData` stores one row per integration and baseline, in the
column layout of ``aas6sim`` ``visibility.py`` (FITS binary table
``VISIBILITY``: TIME, HA, UU, VV, WW, VISREAL, VISIMAG, WEIGHT), extended
by a baseline index and a frequency axis:

=================  ===========================  ==================================
attribute          shape                        meaning
=================  ===========================  ==================================
``time_mjd``       (n_rows,)                    time of the integration (MJD, UTC)
``ha_deg``         (n_rows,)                    hour angle of the phase centre at
                                                the reference station
``u, v, w``        (n_rows,)                    baseline coordinates [wavelengths]
                                                at ``freq_hz[0]``
``baseline``       (n_rows,) int                index into ``baseline_names``
``freq_hz``        (n_chan,)                    channel frequencies
``vis``            (n_rows, n_chan) complex     visibility [Jy]
``weight``         (n_rows, n_chan)             1 / sigma^2 of the real (imaginary) part
=================  ===========================  ==================================

``header`` holds descriptive keywords (telescope, object, RA/Dec,
``SIMULATD`` ...) as a plain dictionary; the FITS input/output follows in
the package ``io``. :meth:`save_npz` / :meth:`load_npz` store the data
losslessly in the meantime.

Version 0.5.0 (2026-10-31): first implementation (design: chapter "Basic
Design", module visdata).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

__all__ = ["VisData"]


@dataclass
class VisData:
    """Visibilities: one row per integration and baseline, n_chan channels."""

    time_mjd: np.ndarray
    ha_deg: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    baseline: np.ndarray
    vis: np.ndarray
    weight: np.ndarray
    freq_hz: np.ndarray
    baseline_names: List[str] = field(default_factory=lambda: ["A-B"])
    header: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.vis = np.asarray(self.vis, dtype=complex)
        if self.vis.ndim == 1:
            self.vis = self.vis[:, None]
        self.weight = np.asarray(self.weight, dtype=float)
        if self.weight.ndim == 1:
            self.weight = self.weight[:, None]
        self.freq_hz = np.atleast_1d(np.asarray(self.freq_hz, dtype=float))
        for name in ("time_mjd", "ha_deg", "u", "v", "w"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        self.baseline = np.asarray(self.baseline, dtype=int)
        n = self.vis.shape[0]
        for name in ("time_mjd", "ha_deg", "u", "v", "w", "baseline"):
            if getattr(self, name).shape != (n,):
                raise ValueError(f"{name}: expected shape ({n},), got {getattr(self, name).shape}")
        if self.weight.shape != self.vis.shape or self.vis.shape[1] != self.freq_hz.size:
            raise ValueError("vis, weight and freq_hz do not match")

    # ------------------------------------------------------------------
    @property
    def n_rows(self) -> int:
        return self.vis.shape[0]

    @property
    def n_chan(self) -> int:
        return self.freq_hz.size

    @property
    def sigma(self) -> np.ndarray:
        """Noise of the real (imaginary) part, 1 / sqrt(weight)."""
        with np.errstate(divide="ignore"):
            return np.where(self.weight > 0, 1.0 / np.sqrt(self.weight), np.inf)

    @property
    def amplitude(self) -> np.ndarray:
        return np.abs(self.vis)

    @property
    def phase_deg(self) -> np.ndarray:
        return np.degrees(np.angle(self.vis))

    def _rows(self, mask: np.ndarray) -> "VisData":
        return replace(self, time_mjd=self.time_mjd[mask], ha_deg=self.ha_deg[mask],
                       u=self.u[mask], v=self.v[mask], w=self.w[mask],
                       baseline=self.baseline[mask], vis=self.vis[mask],
                       weight=self.weight[mask], header=dict(self.header))

    def select_baseline(self, name: Union[str, int]) -> "VisData":
        """Rows of one baseline (by name such as "A-B" or by index)."""
        k = self.baseline_names.index(name) if isinstance(name, str) else int(name)
        return self._rows(self.baseline == k)

    def channel(self, k: int) -> "VisData":
        """One channel (u, v, w are rescaled to its frequency)."""
        s = self.freq_hz[k] / self.freq_hz[0]
        return replace(self, u=self.u * s, v=self.v * s, w=self.w * s,
                       vis=self.vis[:, k:k + 1], weight=self.weight[:, k:k + 1],
                       freq_hz=self.freq_hz[k:k + 1], header=dict(self.header))

    def average_channels(self) -> "VisData":
        """Weighted mean over all channels (one channel at the mean frequency)."""
        wsum = self.weight.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            v = np.where(wsum > 0, (self.vis * self.weight).sum(axis=1) / wsum, 0)
        return replace(self, vis=v[:, None], weight=wsum[:, None],
                       freq_hz=np.array([self.freq_hz.mean()]), header=dict(self.header))

    def concatenate(self, other: "VisData") -> "VisData":
        """Rows of both data sets (same channels and baseline names)."""
        if other.baseline_names != self.baseline_names or not np.allclose(other.freq_hz, self.freq_hz):
            raise ValueError("baseline names or channels differ")
        cat = lambda a, b: np.concatenate([a, b])
        return replace(self, time_mjd=cat(self.time_mjd, other.time_mjd),
                       ha_deg=cat(self.ha_deg, other.ha_deg), u=cat(self.u, other.u),
                       v=cat(self.v, other.v), w=cat(self.w, other.w),
                       baseline=cat(self.baseline, other.baseline),
                       vis=np.concatenate([self.vis, other.vis]),
                       weight=np.concatenate([self.weight, other.weight]),
                       header=dict(self.header))

    # ------------------------------------------------------------------
    def save_npz(self, path: Union[str, Path]) -> Path:
        """Store losslessly as NumPy archive (header as JSON)."""
        path = Path(path)
        np.savez_compressed(path, time_mjd=self.time_mjd, ha_deg=self.ha_deg, u=self.u, v=self.v,
                            w=self.w, baseline=self.baseline, vis=self.vis, weight=self.weight,
                            freq_hz=self.freq_hz,
                            baseline_names=np.array(json.dumps(self.baseline_names)),
                            header=np.array(json.dumps(self.header, default=str)))
        return path if path.suffix == ".npz" else path.with_suffix(path.suffix + ".npz")

    @classmethod
    def load_npz(cls, path: Union[str, Path]) -> "VisData":
        with np.load(path, allow_pickle=False) as d:
            return cls(time_mjd=d["time_mjd"], ha_deg=d["ha_deg"], u=d["u"], v=d["v"], w=d["w"],
                       baseline=d["baseline"], vis=d["vis"], weight=d["weight"],
                       freq_hz=d["freq_hz"],
                       baseline_names=json.loads(str(d["baseline_names"])),
                       header=json.loads(str(d["header"])))

    @classmethod
    def from_measurements(cls, meas: Dict[str, dict], freq_hz: float,
                          start_mjd: float = 0.0, ha0_deg: float = 0.0,
                          header: Optional[Dict[str, Any]] = None) -> "VisData":
        """From the result of ``signal.measure_visibilities`` (notebook 03)."""
        from .observation.geometry import SIDEREAL_DEG_PER_S
        names = list(meas)
        cols = {k: [] for k in ("t", "u", "v", "w", "b", "vis", "wt")}
        for k, name in enumerate(names):
            m = meas[name]
            n = m["vis"].size
            cols["t"].append(m["t_s"]); cols["u"].append(m["uvw"][:, 0])
            cols["v"].append(m["uvw"][:, 1]); cols["w"].append(m["uvw"][:, 2])
            cols["b"].append(np.full(n, k)); cols["vis"].append(m["vis"])
            cols["wt"].append(np.full(n, 1.0 / m["sigma"] ** 2))
        t = np.concatenate(cols["t"])
        return cls(time_mjd=start_mjd + t / 86400.0, ha_deg=ha0_deg + SIDEREAL_DEG_PER_S * t,
                   u=np.concatenate(cols["u"]), v=np.concatenate(cols["v"]),
                   w=np.concatenate(cols["w"]), baseline=np.concatenate(cols["b"]),
                   vis=np.concatenate(cols["vis"]), weight=np.concatenate(cols["wt"]),
                   freq_hz=[freq_hz], baseline_names=names,
                   header=dict(header or {}, SIMULATD=True, ORIGIN="interferosim.signal"))
