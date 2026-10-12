# src/interferosim/sky/models.py -- sky model and JSON format (interferosim.sky.models)
"""Sky model: field centre plus a list of components.

File format
-----------
The JSON format is the one of AAS6 ``json/interfero/sim_sky.json``::

    {
      "_hinweis": "free text",              (or "_note")
      "center": {"label": "W3(OH)", "ra_deg": 36.7661, "dec_deg": 61.8742},
      "components": [
        {"name": "Maser M1", "dx_deg": 0.3, "dy_deg": 0.45, "flux_jy": 70.0,
         "kind": "point", "fwhm_x_deg": 0.0, "fwhm_y_deg": 0.0, "pa_deg": 0.0},
        ...
      ]
    }

Interferosim adds the component kind "disc" (fields ``diameter_deg``,
``limb``, ``limb_power``). Files written by Interferosim keep the AAS6 keys,
so AAS6 can read them as long as only "point" and "gauss" are used.

Version 0.2.0 (2026-10-12): first implementation.
Version 0.2.1 (2026-10-12): demo component names in English.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np

from .components import (GaussianSource, PointSource, SkyComponent,
                         component_from_dict)

__all__ = ["FieldCenter", "SkyModel", "AS", "AM"]

AS = 1.0 / 3600.0          # arc second in degrees
AM = 1.0 / 60.0            # arc minute in degrees


@dataclass
class FieldCenter:
    """Phase centre of the field (J2000, degrees)."""

    label: str = ""
    ra_deg: Optional[float] = None
    dec_deg: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {"label": self.label, "ra_deg": self.ra_deg, "dec_deg": self.dec_deg}

    @classmethod
    def from_dict(cls, d: Optional[Dict[str, Any]]) -> "FieldCenter":
        d = d or {}
        return cls(label=str(d.get("label", "")),
                   ra_deg=None if d.get("ra_deg") is None else float(d["ra_deg"]),
                   dec_deg=None if d.get("dec_deg") is None else float(d["dec_deg"]))


@dataclass
class SkyModel:
    """Field centre plus components; renders images and visibilities."""

    components: List[SkyComponent] = field(default_factory=list)
    center: FieldCenter = field(default_factory=FieldCenter)
    note: str = ""
    name: str = ""

    # ------------------------------------------------------------------
    # building
    def add(self, component: SkyComponent) -> "SkyModel":
        """Append a component; returns the model (allows chaining)."""
        if not isinstance(component, SkyComponent):
            raise TypeError("component must be a SkyComponent")
        self.components.append(component)
        return self

    def __len__(self) -> int:
        return len(self.components)

    def __iter__(self):
        return iter(self.components)

    def __getitem__(self, name: str) -> SkyComponent:
        for c in self.components:
            if c.name == name:
                return c
        raise KeyError(name)

    @property
    def total_flux_jy(self) -> float:
        return float(sum(c.flux_jy for c in self.components))

    # ------------------------------------------------------------------
    # geometry of the field
    def bounding_box_deg(self) -> Tuple[float, float, float, float]:
        """(x_min, x_max, y_min, y_max) of all components incl. their extent."""
        if not self.components:
            return (-0.5, 0.5, -0.5, 0.5)
        h = [0.5 * c.extent_deg for c in self.components]
        xs = [c.dx_deg for c in self.components]
        ys = [c.dy_deg for c in self.components]
        return (min(x - e for x, e in zip(xs, h)), max(x + e for x, e in zip(xs, h)),
                min(y - e for y, e in zip(ys, h)), max(y + e for y, e in zip(ys, h)))

    def map_center_deg(self) -> Tuple[float, float]:
        """Centre of the bounding box (default centre of a map)."""
        x0, x1, y0, y1 = self.bounding_box_deg()
        return 0.5 * (x0 + x1), 0.5 * (y0 + y1)

    def extent_deg(self, margin: float = 1.4) -> float:
        """Full width of a square map around :meth:`map_center_deg` that
        contains all components."""
        if not self.components:
            return 1.0
        x0, x1, y0, y1 = self.bounding_box_deg()
        width = max(x1 - x0, y1 - y0)
        if width <= 0:                    # only point sources at one position
            width = 60.0 / 3600.0
        return float(width * margin)

    def grid(self, fov_deg: Optional[float] = None, pix_deg: Optional[float] = None,
             npix: int = 512, center_deg: Optional[Tuple[float, float]] = None):
        """Pixel axes (x_deg, y_deg) and pixel size of a square map.

        ``fov_deg`` defaults to :meth:`extent_deg`; ``pix_deg`` to fov/npix.
        ``center_deg`` is the map centre as offset from the field centre;
        default :meth:`map_center_deg`, so that all components are in the
        map even if they are not grouped around the field centre.
        """
        if center_deg is None:
            center_deg = self.map_center_deg()
        fov = float(fov_deg) if fov_deg else self.extent_deg()
        pix = float(pix_deg) if pix_deg else fov / npix
        n = max(int(round(fov / pix)), 1)
        ax = (np.arange(n) - n // 2) * pix
        return ax + center_deg[0], ax + center_deg[1], pix

    # ------------------------------------------------------------------
    # products
    def render(self, x_deg: np.ndarray, y_deg: np.ndarray, pix_deg: float,
               dtype=np.float32) -> np.ndarray:
        """Image in Jy per pixel on the axes x_deg (east), y_deg (north)."""
        img = np.zeros((np.size(y_deg), np.size(x_deg)), dtype=dtype)
        for c in self.components:
            c.render(x_deg, y_deg, pix_deg, img)
        return img

    def render_map(self, fov_deg: Optional[float] = None, pix_deg: Optional[float] = None,
                   npix: int = 512, center_deg: Optional[Tuple[float, float]] = None):
        """Convenience: returns (image, x_deg, y_deg, pix_deg)."""
        x, y, pix = self.grid(fov_deg, pix_deg, npix, center_deg)
        return self.render(x, y, pix), x, y, pix

    def visibility(self, u, v) -> np.ndarray:
        """Sum of the analytic component visibilities in Jy."""
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        V = np.zeros(np.broadcast(u, v).shape, dtype=complex)
        for c in self.components:
            V = V + c.visibility(u, v)
        return V

    # ------------------------------------------------------------------
    # file format (AAS6 sim_sky.json)
    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {}
        if self.note:
            d["_note"] = self.note
        if self.name:
            d["name"] = self.name
        d["center"] = self.center.to_dict()
        d["components"] = [c.to_dict() for c in self.components]
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SkyModel":
        comps = [component_from_dict(c) for c in d.get("components", [])]
        note = str(d.get("_note", d.get("_hinweis", "")))
        return cls(components=comps, center=FieldCenter.from_dict(d.get("center")),
                   note=note, name=str(d.get("name", "")))

    def to_json(self, path: Union[str, Path], indent: int = 2) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=indent, ensure_ascii=False),
                        encoding="utf-8")
        return path

    @classmethod
    def from_json(cls, path: Union[str, Path]) -> "SkyModel":
        path = Path(path)
        model = cls.from_dict(json.loads(path.read_text(encoding="utf-8")))
        if not model.name:
            model.name = path.stem
        return model

    # ------------------------------------------------------------------
    @classmethod
    def demo(cls) -> "SkyModel":
        """Demo field of aas6sim (``model.demo_sky()``), centred on W3(OH)
        as in AAS6 ``json/interfero/sim_sky.json``.

        Positions and fluxes are identical to aas6sim; since version 0.2.1 the
        component names are English (project language). Names in aas6sim /
        AAS6: HI-Sockel = HI base, Gal. Riegel = Galactic bar, HII-Region A/B =
        HII region A/B, SNR-Schale C = SNR shell C, Kompaktquelle K1 =
        Compact source K1; the masers keep their names.
        """
        comps: Iterable[SkyComponent] = [
            GaussianSource("HI base", 0.0, 0.0, 1500.0, 12.0, 12.0),
            GaussianSource("Galactic bar", 0.0, 0.0, 4000.0, 8.0, 1.6, pa_deg=25.0),
            GaussianSource("HII region A", -2.0, 1.2, 900.0, 0.70, 0.70),
            GaussianSource("HII region B", 2.6, -1.6, 550.0, 0.40, 0.40),
            GaussianSource("SNR shell C", 1.2, 2.4, 350.0, 1.10, 0.85, pa_deg=-40.0),
            GaussianSource("Compact source K1", -1.4, -0.6, 260.0, 2.1 * AM, 2.1 * AM),
            PointSource("Maser M1", 0.30, 0.45, 70.0),
            PointSource("Maser M2", 0.30 + 5.2 * AS, 0.45 - 1.5 * AS, 45.0),
            PointSource("Maser M3", 0.30 + 14.0 * AS, 0.45 + 9.0 * AS, 20.0),
        ]
        return cls(components=list(comps),
                   center=FieldCenter("W3(OH)", 36.7661, 61.8742),
                   note="Demo field of aas6sim (model.demo_sky()), field centre W3(OH).",
                   name="aas6_demo")
