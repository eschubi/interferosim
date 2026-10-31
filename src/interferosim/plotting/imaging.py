# src/interferosim/plotting/imaging.py -- image plots (interferosim.plotting.imaging)
"""Plots of images (dirty image, dirty beam, CLEAN results).

* :func:`plot_image` -- map in arcsec offsets from the image centre, east
  to the left, optional zoom, source marks and colour bar.
* :func:`sky_marks` -- marks (l, m, name) of the components of a sky model
  relative to an image centre, for :func:`plot_image`.
* :func:`plot_clean_history` -- convergence of a CLEAN run: residual peak
  and cumulative component flux against the iteration (since 0.8.0).

Version 0.8.0 (2026-10-31): plot_clean_history.
Version 0.7.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ..imaging import Image
from ..sky import SkyModel

__all__ = ["plot_image", "sky_marks", "plot_clean_history"]


def sky_marks(sky: SkyModel, center_deg: Tuple[float, float] = (0.0, 0.0),
              min_flux_jy: float = 0.0) -> List[Tuple[float, float, str]]:
    """(l, m, name) in arcsec of all components above ``min_flux_jy``."""
    return [((c.dx_deg - center_deg[0]) * 3600.0, (c.dy_deg - center_deg[1]) * 3600.0, c.name)
            for c in sky if c.flux_jy >= min_flux_jy]


def plot_image(img: Image, *, ax: Optional[plt.Axes] = None, zoom_arcsec: Optional[float] = None,
               marks: Optional[Iterable[Tuple[float, float, str]]] = None, cmap: str = "inferno",
               vmin: Optional[float] = None, vmax: Optional[float] = None, colorbar: bool = True,
               title: Optional[str] = None, labels: bool = False):
    """Show ``img``; ``zoom_arcsec`` is the half width of the shown region.
    Returns (fig, ax)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.2, 5.2))
    else:
        fig = ax.figure
    im = ax.imshow(img.data, origin="lower", extent=img.extent_arcsec(), cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    if zoom_arcsec:
        ax.set_xlim(zoom_arcsec, -zoom_arcsec); ax.set_ylim(-zoom_arcsec, zoom_arcsec)
    else:
        e = img.extent_arcsec()
        ax.set_xlim(e[1], e[0])                                   # east to the left
    for l, m, name in (marks or []):
        ax.plot(l, m, "o", ms=11, mfc="none", mec="cyan", mew=1.2)
        if labels:
            ax.annotate(name, (l, m), xytext=(7, 7), textcoords="offset points", color="cyan", fontsize=8)
    ax.set_xlabel(r"$\Delta\alpha$ (east)  [arcsec]"); ax.set_ylabel(r"$\Delta\delta$ (north)  [arcsec]")
    if colorbar:
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(img.unit or "relative")
    ax.set_title(title or img.name, fontsize=10)
    return fig, ax


def plot_clean_history(result, *, ax: Optional[plt.Axes] = None, label: Optional[str] = None,
                       title: Optional[str] = None):
    """Residual peak (left axis, log) and cumulative CLEAN flux (right axis)
    against the iteration of a :class:`~interferosim.deconvolution.DeconvResult`.
    Returns (fig, ax, ax_right)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
    else:
        fig = ax.figure
    h = result.history
    it = np.arange(1, len(h.get("peak", [])) + 1)
    ax.semilogy(it, h["peak"], color="C0", label=label or "residual peak")
    thr = result.meta.get("threshold_jy", 0.0)
    if thr:
        ax.axhline(thr, color="C0", ls=":", lw=1.0)
    ax.set_xlabel("iteration"); ax.set_ylabel("residual peak  [Jy/beam]", color="C0")
    ax.grid(alpha=0.3, which="both")
    ax2 = ax.twinx()
    ax2.plot(it, h["flux"], color="C3")
    ax2.set_ylabel("flux in CLEAN components  [Jy]", color="C3")
    ax.set_title(title or f"CLEAN convergence, {result.n_iter} iterations", fontsize=10)
    return fig, ax, ax2
