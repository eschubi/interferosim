# src/interferosim/plotting/observation.py -- observation plots (interferosim.plotting.observation)
"""Plots of observations.

* :func:`plot_uv_coverage` -- uv-tracks of all baselines (and their
  conjugates -V(-u, -v)), visible part solid, below the elevation limit grey.
* :func:`plot_elevation` -- elevation of the source over the hour angle at
  every station, with the elevation limit.
* :func:`plot_visdata` -- amplitude and phase of a VisData set over the
  hour angle, per baseline.

Version 0.5.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..observation import Observation, UVCoverage
from ..visdata import VisData

__all__ = ["plot_uv_coverage", "plot_elevation", "plot_visdata"]


def plot_uv_coverage(cov: UVCoverage, *, ax: Optional[plt.Axes] = None, unit: str = "auto",
                     title: Optional[str] = None):
    """uv-coverage; ``unit`` "lambda", "klambda" or "m" (auto: kilolambda
    above 10^4 wavelengths). Returns (fig, ax)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.2, 6.0))
    else:
        fig = ax.figure
    qmax = max(float(np.max(np.hypot(cov.u, cov.v))), 1e-9)
    if unit == "auto":
        unit = "klambda" if qmax > 1e4 else "lambda"
    f = {"lambda": 1.0, "klambda": 1e-3, "m": cov.wavelength_m}[unit]
    lab = {"lambda": r"[$\lambda$]", "klambda": r"[k$\lambda$]", "m": "[m]"}[unit]
    for k, name in enumerate(cov.baseline_names):
        u, v = cov.u[k] * f, cov.v[k] * f
        line, = ax.plot(np.where(cov.mask, u, np.nan), np.where(cov.mask, v, np.nan), lw=1.6, label=name)
        c = line.get_color()
        ax.plot(np.where(cov.mask, -u, np.nan), np.where(cov.mask, -v, np.nan), lw=1.6, color=c, alpha=0.5)
        ax.plot(np.where(~cov.mask, u, np.nan), np.where(~cov.mask, v, np.nan), lw=0.6, color="grey")
        ax.plot(np.where(~cov.mask, -u, np.nan), np.where(~cov.mask, -v, np.nan), lw=0.6, color="grey")
    ax.set_xlabel(f"u (east)  {lab}"); ax.set_ylabel(f"v (north)  {lab}")
    ax.set_aspect("equal", adjustable="datalim"); ax.grid(alpha=0.3)
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    if len(cov.baseline_names) <= 8:
        ax.legend(fontsize=8, loc="upper right")
    ax.set_title(title or f"uv-coverage, resolution ~{cov.resolution_arcsec():.3g}\"", fontsize=10)
    return fig, ax


def plot_elevation(obs: Observation, *, ax: Optional[plt.Axes] = None, ha_range_h=(-12.0, 12.0)):
    """Elevation over the hour angle at every station. Returns (fig, ax)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.0, 4.0))
    else:
        fig = ax.figure
    ha = np.linspace(ha_range_h[0] * 15, ha_range_h[1] * 15, 721)
    el = obs.elevation_deg(ha)
    for k, st in enumerate(obs.array.stations):
        ax.plot(ha / 15.0, el[k], lw=1.2, label=f"station {st.name}")
    ax.axhline(obs.min_elev_deg, color="C3", ls="--", lw=1, label=f"limit {obs.min_elev_deg:g} deg")
    ax.axvspan(*obs.ha_range_h, color="C2", alpha=0.08, label="observed range")
    ax.set_xlabel("hour angle  [h]"); ax.set_ylabel("elevation  [deg]")
    ax.set_ylim(-10, 90); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"{obs.name or 'source'} (Dec {obs.dec_deg:+.2f} deg) at {obs.array.site.name}",
                 fontsize=10)
    return fig, ax


def plot_visdata(vd: VisData, *, channel: int = 0, model: Optional[VisData] = None,
                 figsize=(12.0, 6.0)):
    """Amplitude and phase over the hour angle per baseline; optional model
    (e.g. without noise) as line. Returns (fig, axes)."""
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    for k, name in enumerate(vd.baseline_names):
        sel = vd.baseline == k
        h = vd.ha_deg[sel] / 15.0
        p = axes[0].plot(h, np.abs(vd.vis[sel, channel]), ".", ms=3, label=name)
        c = p[0].get_color()
        axes[1].plot(h, np.degrees(np.angle(vd.vis[sel, channel])), ".", ms=3, color=c)
        if model is not None:
            ms = model.baseline == k
            axes[0].plot(model.ha_deg[ms] / 15.0, np.abs(model.vis[ms, channel]), "-", lw=1, color=c)
            axes[1].plot(model.ha_deg[ms] / 15.0, np.degrees(np.angle(model.vis[ms, channel])),
                         "-", lw=0.8, color=c, alpha=0.7)
    axes[0].set_ylabel("|V|  [Jy]"); axes[1].set_ylabel("phase  [deg]"); axes[1].set_ylim(-190, 190)
    axes[1].set_xlabel("hour angle  [h]")
    for ax in axes:
        ax.grid(alpha=0.3)
    if len(vd.baseline_names) <= 8:
        axes[0].legend(fontsize=8)
    fig.suptitle(f"{vd.header.get('OBJECT', '')}  ({vd.freq_hz[channel] / 1e6:.3f} MHz)", fontsize=10)
    fig.tight_layout()
    return fig, axes
