# src/interferosim/sky/components.py -- sky components (interferosim.sky.components)
"""Sky components of the Interferosim sky model.

Every component knows its position relative to the field centre, its flux
density, how it is rendered onto a pixel grid and its analytic visibility.

Origin
------
Point and Gaussian components and their rendering/visibility are taken
over from aas6sim 1.3.1 (``model.Component``, ``model.render``,
``model.visibility``) and split into one class per source type. The disc
component (uniform or limb-brightened, e.g. the Sun at 21 cm) comes from
notebook B (Appendix B, ``visibility_profile`` and ``PROFILES``).

Conventions
-----------
* ``dx_deg``/``dy_deg``: offset from the field centre in degrees, east and
  north (same keys as ``json/interfero/sim_sky.json`` in AAS6). In the
  small-field approximation they are the direction cosines (l, m).
* ``flux_jy``: total flux density in Jy.
* Rendered images are in Jy per pixel; ``u``, ``v`` are in wavelengths.
* ``pa_deg``: rotation of the Gaussian major axis (``fwhm_x_deg``) from the
  +dx axis towards +dy -- the convention of aas6sim, kept unchanged.

Version 0.2.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, fields
from typing import Any, ClassVar, Dict, Type

import numpy as np
from scipy.special import j0, j1

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
# exp(-GAUSS_UV * (theta*q)^2) is the visibility of a Gaussian of FWHM theta
GAUSS_UV = np.pi**2 / (4.0 * np.log(2.0))

__all__ = [
    "SkyComponent",
    "PointSource",
    "GaussianSource",
    "DiscSource",
    "component_from_dict",
    "COMPONENT_KINDS",
]


def _nearest_index(axis: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(axis - value)))


@dataclass
class SkyComponent(ABC):
    """Abstract base class of all sky components."""

    name: str
    dx_deg: float
    dy_deg: float
    flux_jy: float

    kind: ClassVar[str] = ""

    # ------------------------------------------------------------------
    @property
    def is_point(self) -> bool:
        """True if the component has no extent."""
        return False

    @property
    def extent_deg(self) -> float:
        """Approximate full size (for automatic field-of-view selection)."""
        return 0.0

    # ------------------------------------------------------------------
    @abstractmethod
    def render(self, x_deg: np.ndarray, y_deg: np.ndarray, pix_deg: float,
               out: np.ndarray) -> None:
        """Add the component to ``out`` (shape (ny, nx), Jy per pixel)."""

    @abstractmethod
    def visibility(self, u, v) -> np.ndarray:
        """Analytic visibility in Jy for spatial frequencies u, v [lambda]."""

    # ------------------------------------------------------------------
    def _phase(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        l0, m0 = np.radians(self.dx_deg), np.radians(self.dy_deg)
        return np.exp(-2j * np.pi * (u * l0 + v * m0))

    def to_dict(self) -> Dict[str, Any]:
        """Dictionary in the sim_sky.json component format."""
        d = {"name": self.name, "dx_deg": self.dx_deg, "dy_deg": self.dy_deg,
             "flux_jy": self.flux_jy, "kind": self.kind}
        for k, v in asdict(self).items():
            if k not in d:
                d[k] = v
        return d


# ======================================================================
@dataclass
class PointSource(SkyComponent):
    """Unresolved source (kind "point")."""

    kind: ClassVar[str] = "point"

    @property
    def is_point(self) -> bool:
        return True

    def render(self, x_deg, y_deg, pix_deg, out) -> None:
        i = _nearest_index(np.asarray(x_deg, float), self.dx_deg)
        j = _nearest_index(np.asarray(y_deg, float), self.dy_deg)
        out[j, i] += self.flux_jy

    def visibility(self, u, v) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        return self.flux_jy * self._phase(u, v)

    def to_dict(self) -> Dict[str, Any]:
        # sim_sky.json stores the Gaussian fields with 0 for point sources
        d = super().to_dict()
        d.update({"fwhm_x_deg": 0.0, "fwhm_y_deg": 0.0, "pa_deg": 0.0})
        return d


# ======================================================================
@dataclass
class GaussianSource(SkyComponent):
    """Elliptical Gaussian source (kind "gauss").

    ``fwhm_y_deg`` <= 0 means circular (same as ``fwhm_x_deg``).
    """

    fwhm_x_deg: float = 0.0
    fwhm_y_deg: float = 0.0
    pa_deg: float = 0.0

    kind: ClassVar[str] = "gauss"

    @property
    def fwhm_pair_deg(self):
        fx = self.fwhm_x_deg
        return fx, (self.fwhm_y_deg if self.fwhm_y_deg > 0 else fx)

    @property
    def is_point(self) -> bool:
        return max(self.fwhm_x_deg, self.fwhm_y_deg) <= 0.0

    @property
    def extent_deg(self) -> float:
        return 2.0 * max(self.fwhm_pair_deg)

    def render(self, x_deg, y_deg, pix_deg, out, chunk: int = 128) -> None:
        # Algorithm of aas6sim 1.3.1 model.render: components smaller than
        # two pixels are put into a single pixel; otherwise the Gaussian is
        # evaluated in blocks of `chunk` rows to keep memory use small.
        x = np.asarray(x_deg, dtype=float)
        y = np.asarray(y_deg, dtype=float)
        if self.is_point or max(self.fwhm_pair_deg) < 2 * pix_deg:
            out[_nearest_index(y, self.dy_deg), _nearest_index(x, self.dx_deg)] += self.flux_jy
            return
        fx, fy = self.fwhm_pair_deg
        sx, sy = fx * FWHM_TO_SIGMA, fy * FWHM_TO_SIGMA
        pa = np.radians(self.pa_deg)
        cp, sp = np.cos(pa), np.sin(pa)
        amp = self.flux_jy / (2 * np.pi * sx * sy) * pix_deg**2
        dx = x - self.dx_deg
        for j0_ in range(0, y.size, chunk):
            dy = (y[j0_:j0_ + chunk] - self.dy_deg)[:, None]
            xr = (dx[None, :] * cp + dy * sp) / sx
            yr = (dx[None, :] * (-sp) + dy * cp) / sy
            out[j0_:j0_ + chunk] += (amp * np.exp(-0.5 * (xr * xr + yr * yr))).astype(
                out.dtype, copy=False)

    def visibility(self, u, v) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        if self.is_point:
            return self.flux_jy * self._phase(u, v)
        fx, fy = self.fwhm_pair_deg
        tx, ty = np.radians(fx), np.radians(fy)
        pa = np.radians(self.pa_deg)
        ur = u * np.cos(pa) + v * np.sin(pa)
        vr = -u * np.sin(pa) + v * np.cos(pa)
        env = np.exp(-GAUSS_UV * ((tx * ur) ** 2 + (ty * vr) ** 2))
        return self.flux_jy * env * self._phase(u, v)


# ======================================================================
@dataclass
class DiscSource(SkyComponent):
    """Circular disc (kind "disc"), uniform or limb-brightened.

    Radial brightness profile (notebook B)::

        I(rho) = 1 + limb * rho**limb_power,   rho = r / R  (0 <= rho <= 1)

    ``limb = 0`` gives the uniform disc with the analytic visibility
    2 J1(x)/x; otherwise the visibility is the numerical Hankel transform.
    """

    diameter_deg: float = 0.0
    limb: float = 0.0
    limb_power: float = 8.0

    kind: ClassVar[str] = "disc"

    @property
    def is_point(self) -> bool:
        return self.diameter_deg <= 0.0

    @property
    def extent_deg(self) -> float:
        return self.diameter_deg

    def profile(self, rho: np.ndarray) -> np.ndarray:
        """Relative brightness for rho = r/R (zero outside the disc)."""
        rho = np.asarray(rho, dtype=float)
        p = 1.0 + self.limb * np.clip(rho, 0.0, 1.0) ** self.limb_power
        return np.where(rho <= 1.0, p, 0.0)

    def render(self, x_deg, y_deg, pix_deg, out) -> None:
        x = np.asarray(x_deg, dtype=float)
        y = np.asarray(y_deg, dtype=float)
        r_deg = 0.5 * self.diameter_deg
        if self.is_point or self.diameter_deg < 2 * pix_deg:
            out[_nearest_index(y, self.dy_deg), _nearest_index(x, self.dx_deg)] += self.flux_jy
            return
        # only the bounding box of the disc is evaluated
        ix = np.nonzero(np.abs(x - self.dx_deg) <= r_deg + pix_deg)[0]
        iy = np.nonzero(np.abs(y - self.dy_deg) <= r_deg + pix_deg)[0]
        if ix.size == 0 or iy.size == 0:
            return
        X, Y = np.meshgrid(x[ix] - self.dx_deg, y[iy] - self.dy_deg)
        img = self.profile(np.hypot(X, Y) / r_deg)
        total = img.sum()
        if total <= 0:
            return
        out[iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1] += (img * (self.flux_jy / total)).astype(
            out.dtype, copy=False)

    def visibility(self, u, v, n_r: int = 2000) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        q = np.hypot(u, v)
        if self.is_point:
            return self.flux_jy * self._phase(u, v)
        theta = np.radians(self.diameter_deg)
        if self.limb == 0.0:
            xq = np.pi * theta * q
            with np.errstate(invalid="ignore", divide="ignore"):
                env = np.where(xq < 1e-9, 1.0, 2.0 * j1(xq) / np.where(xq == 0, 1, xq))
        else:
            # Hankel transform, normalised to V(0) = 1 (notebook B)
            R = 0.5 * theta
            r = np.linspace(0.0, R, n_r)
            w = self.profile(r / R) * r
            kern = j0(2.0 * np.pi * np.outer(q.ravel(), r))
            num = np.trapezoid(kern * w, r, axis=-1) if hasattr(np, "trapezoid") \
                else np.trapz(kern * w, r, axis=-1)
            den = np.trapezoid(w, r) if hasattr(np, "trapezoid") else np.trapz(w, r)
            env = (num / den).reshape(q.shape)
        return self.flux_jy * env * self._phase(u, v)


# ======================================================================
COMPONENT_KINDS: Dict[str, Type[SkyComponent]] = {
    PointSource.kind: PointSource,
    GaussianSource.kind: GaussianSource,
    DiscSource.kind: DiscSource,
}


def component_from_dict(d: Dict[str, Any]) -> SkyComponent:
    """Create a component from a sim_sky.json component entry.

    Compatible with AAS6: "kind" defaults to "gauss" (as in aas6sim), and a
    "gauss" entry without extent is treated as a point source there; here it
    stays a GaussianSource whose ``is_point`` is True, which behaves the same.
    """
    kind = str(d.get("kind", "gauss")).lower()
    if kind not in COMPONENT_KINDS:
        raise ValueError(f"unknown component kind {kind!r} "
                         f"(known: {', '.join(COMPONENT_KINDS)})")
    cls = COMPONENT_KINDS[kind]
    allowed = {f.name for f in fields(cls)}          # ClassVar 'kind' is not a field
    kwargs = {k: v for k, v in d.items() if k in allowed}
    missing = {"name", "dx_deg", "dy_deg", "flux_jy"} - set(kwargs)
    if missing:
        raise ValueError(f"component {d.get('name', '?')!r}: missing {sorted(missing)}")
    return cls(**kwargs)
