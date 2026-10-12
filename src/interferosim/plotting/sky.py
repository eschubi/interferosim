# src/interferosim/plotting/sky.py -- sky maps (interferosim.plotting.sky)
"""Maps of a :class:`~interferosim.sky.SkyModel`.

* :func:`plot_sky_map` -- one map: rendered brightness (Jy per pixel) on a
  linear, logarithmic or asinh scale, point sources and compact components
  marked with circles whose area scales with the flux density (as in the
  AAS6 star-map overlay), optional labels.
* :func:`plot_sky_overview` -- full field plus a zoom on a chosen region
  (e.g. the maser group, which is far below the pixel size of the full map).

Axes: offset from the field centre in degrees; east is to the LEFT
(astronomical convention, as in the notebooks).

Version 0.2.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

from ..sky import SkyModel

__all__ = ["plot_sky_map", "plot_sky_overview"]

_SCALES = ("log", "linear", "asinh")


def _norm(img: np.ndarray, scale: str, dynamic_range: float):
    vmax = float(np.nanmax(img)) if img.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    if scale == "linear":
        return colors.Normalize(vmin=0.0, vmax=vmax)
    if scale == "log":
        return colors.LogNorm(vmin=vmax / dynamic_range, vmax=vmax, clip=True)
    if scale == "asinh":
        if hasattr(colors, "AsinhNorm"):
            return colors.AsinhNorm(linear_width=vmax / dynamic_range, vmin=0.0, vmax=vmax)
        return colors.PowerNorm(gamma=0.3, vmin=0.0, vmax=vmax)
    raise ValueError(f"scale must be one of {_SCALES}")


def plot_sky_map(sky: SkyModel, *, fov_deg: Optional[float] = None,
                 pix_deg: Optional[float] = None, npix: int = 512,
                 center_deg: Optional[Tuple[float, float]] = None,
                 scale: str = "log", dynamic_range: float = 1e4,
                 cmap: str = "inferno", mark: bool = True, labels: bool = True,
                 mark_max_size_deg: Optional[float] = None,
                 ax: Optional[plt.Axes] = None, title: Optional[str] = None,
                 colorbar: bool = True):
    """Render ``sky`` and show it as a map. Returns (fig, ax, image).

    Parameters
    ----------
    fov_deg, pix_deg, npix, center_deg
        Map geometry, see :meth:`SkyModel.grid` (default: all components).
    scale
        "log" (default), "linear" or "asinh".
    dynamic_range
        Ratio max/min of the colour scale for "log" and "asinh".
    mark
        Mark components whose size is below ``mark_max_size_deg``
        (default: 4 pixels) -- they would otherwise be invisible.
    labels
        Write the component names next to the marks.
    """
    if scale not in _SCALES:
        raise ValueError(f"scale must be one of {_SCALES}")
    img, x, y, pix = sky.render_map(fov_deg, pix_deg, npix, center_deg)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.2, 6.2))
    else:
        fig = ax.figure

    floor = float(np.nanmax(img)) / dynamic_range if img.size and np.nanmax(img) > 0 else 1e-30
    shown = np.where(img > 0, img, floor * 1e-3) if scale == "log" else img
    extent = (x[0] - pix / 2, x[-1] + pix / 2, y[0] - pix / 2, y[-1] + pix / 2)
    # axis unit follows the field size: deg, arcmin or arcsec
    fov = extent[1] - extent[0]
    unit, f = ("deg", 1.0) if fov >= 2.0 else (("arcmin", 60.0) if fov >= 0.1
                                                else ("arcsec", 3600.0))
    ext_u = tuple(e * f for e in extent)
    im = ax.imshow(shown, origin="lower", extent=ext_u, cmap=cmap,
                   norm=_norm(img, scale, dynamic_range), interpolation="nearest")
    ax.set_xlim(ext_u[1], ext_u[0])                   # east to the left
    ax.set_ylim(ext_u[2], ext_u[3])
    ax.set_xlabel(rf"$\Delta\alpha$ (east)  [{unit}]")
    ax.set_ylabel(rf"$\Delta\delta$ (north)  [{unit}]")
    ax.set_aspect("equal")

    if mark:
        limit = mark_max_size_deg if mark_max_size_deg is not None else 4 * pix
        small = [c for c in sky.components
                 if c.extent_deg < limit
                 and extent[0] <= c.dx_deg <= extent[1]
                 and extent[2] <= c.dy_deg <= extent[3]]
        if small:
            fmax = max(c.flux_jy for c in small)
            # label offsets alternate so that close groups stay readable
            offsets = [(8, 6), (8, -16), (-60, 6), (-60, -16)]
            for k, c in enumerate(sorted(small, key=lambda c: -c.flux_jy)):
                s = 40.0 + 260.0 * c.flux_jy / fmax          # marker area ~ flux
                ax.scatter([c.dx_deg * f], [c.dy_deg * f], s=s, facecolors="none",
                           edgecolors="cyan", linewidths=1.3, zorder=3)
                if c.is_point:                               # exact position
                    ax.plot([c.dx_deg * f], [c.dy_deg * f], "+", color="white", ms=5, zorder=4)
                if labels:
                    ax.annotate(f"{c.name}\n{c.flux_jy:g} Jy", (c.dx_deg * f, c.dy_deg * f),
                                xytext=offsets[k % len(offsets)], textcoords="offset points",
                                color="cyan", fontsize=7.5, zorder=4)

    if colorbar:
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label("brightness  [Jy / pixel]")
    if title is None:
        c = sky.center
        pos = (f"  (RA {c.ra_deg:.4f}°, Dec {c.dec_deg:+.4f}°)"
               if c.ra_deg is not None and c.dec_deg is not None else "")
        head = sky.name or "sky model"
        title = (f"{head}: centre {c.label or '-'}{pos}\n"
                 f"{len(sky)} components, total {sky.total_flux_jy:,.0f} Jy, "
                 f"pixel {pix * 3600:.3g}\"")
    ax.set_title(title, fontsize=9.5)
    return fig, ax, img


def plot_sky_overview(sky: SkyModel, *, zoom_center_deg: Optional[Tuple[float, float]] = None,
                      zoom_fov_deg: Optional[float] = None, npix: int = 512,
                      scale: str = "log", dynamic_range: float = 1e4,
                      figsize: Sequence[float] = (14.0, 6.0)):
    """Full field (left) and zoom (right). Returns (fig, axes).

    Without ``zoom_center_deg`` the zoom is centred on the brightest point
    source (or the brightest component); ``zoom_fov_deg`` defaults to
    1/20 of the full field.
    """
    full = sky.extent_deg()
    if zoom_center_deg is None:
        cand = [c for c in sky.components if c.is_point] or list(sky.components)
        best = max(cand, key=lambda c: c.flux_jy) if cand else None
        zoom_center_deg = (best.dx_deg, best.dy_deg) if best else (0.0, 0.0)
    if zoom_fov_deg is None:
        # if further point sources lie close to the zoom centre (e.g. a maser
        # group), size the zoom to show all of them; else 1/20 of the field
        near = [np.hypot(c.dx_deg - zoom_center_deg[0], c.dy_deg - zoom_center_deg[1])
                for c in sky.components if c.is_point]
        near = [d for d in near if 0 < d < full / 20.0]
        zoom_fov_deg = 6.0 * max(near) if near else full / 20.0
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    plot_sky_map(sky, fov_deg=full, npix=npix, scale=scale,
                 dynamic_range=dynamic_range, ax=axes[0], labels=False)
    zc = zoom_center_deg
    h = zoom_fov_deg / 2
    lab = axes[0].get_xlabel()                    # unit of the full map
    f0 = 1.0 if "[deg]" in lab else (60.0 if "[arcmin]" in lab else 3600.0)
    axes[0].add_patch(plt.Rectangle(((zc[0] - h) * f0, (zc[1] - h) * f0),
                                    zoom_fov_deg * f0, zoom_fov_deg * f0,
                                    fill=False, edgecolor="white", lw=1.0, ls="--"))
    plot_sky_map(sky, fov_deg=zoom_fov_deg, npix=npix, center_deg=zc, scale=scale,
                 dynamic_range=dynamic_range, ax=axes[1],
                 title=f"zoom: {zoom_fov_deg * 60:.3g}' around "
                       f"({zc[0]:+.4f}°, {zc[1]:+.4f}°)")
    fig.tight_layout()
    return fig, axes
