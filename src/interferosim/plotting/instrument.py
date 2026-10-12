# src/interferosim/plotting/instrument.py -- instrument plots (interferosim.plotting.instrument)
"""Plots of antennas and arrays.

* :func:`plot_beams` -- normalised power patterns of several antennas,
  linear and in dB, with the HPBW marked.
* :func:`plot_array` -- station layout in the local frame (east, north)
  with all baselines and their lengths.

Version 0.3.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from ..constants import HI_21CM
from ..instrument import Antenna, Array

__all__ = ["plot_beams", "plot_array"]


def plot_beams(antennas: Sequence[Antenna], f_hz: float = HI_21CM, *,
               theta_max_deg: Optional[float] = None, figsize=(12.0, 4.5)):
    """Power patterns of ``antennas`` at ``f_hz``. Returns (fig, axes).

    Left: linear, right: in dB (down to -40 dB). ``theta_max_deg`` defaults
    to 1.5 times the largest HPBW.
    """
    antennas = list(antennas)
    if theta_max_deg is None:
        theta_max_deg = 1.5 * max(a.hpbw_deg(f_hz) for a in antennas)
    th = np.linspace(-theta_max_deg, theta_max_deg, 2001)
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for a in antennas:
        p = a.beam(np.abs(th), f_hz)
        h = a.hpbw_deg(f_hz)
        lab = f"{a.name}  ({a.gain_dbi(f_hz):.1f} dBi, HPBW {h:.2f}°)"
        line, = axes[0].plot(th, p, label=lab)
        axes[1].plot(th, 10 * np.log10(np.maximum(p, 1e-12)), color=line.get_color())
        for ax, y in ((axes[0], 0.5), (axes[1], -3.0)):
            ax.plot([-h / 2, h / 2], [y, y], "|", color=line.get_color(), ms=10)
    axes[0].set_ylabel("normalised power")
    axes[1].set_ylabel("power  [dB]")
    axes[1].set_ylim(-40, 3)
    for ax in axes:
        ax.set_xlabel("angle from the axis  [deg]")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle(f"primary beams at {f_hz / 1e6:.3f} MHz", fontsize=10)
    fig.tight_layout()
    return fig, axes


def plot_array(array: Array, *, ax: Optional[plt.Axes] = None, labels: bool = True):
    """Station layout (east to the right, north up) with all baselines.
    Returns (fig, ax)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
    else:
        fig = ax.figure
    enu = array.station_enu()
    span = float(np.ptp(enu[:, :2], axis=0).max()) if len(enu) > 1 else 1.0
    unit, f = ("km", 1e-3) if span >= 2000 else ("m", 1.0)
    for (i, j), L in zip(array.baseline_pairs(), array.baseline_lengths_m()):
        ax.plot(enu[[i, j], 0] * f, enu[[i, j], 1] * f, "-", color="grey", lw=1, zorder=1)
        if labels:
            xm, ym = enu[[i, j], 0].mean() * f, enu[[i, j], 1].mean() * f
            ax.annotate(f"{L * f:.4g} {unit}", (xm, ym), textcoords="offset points",
                        xytext=(0, 6), ha="center", fontsize=8, color="grey")
    ax.scatter(enu[:, 0] * f, enu[:, 1] * f, s=60, color="tab:blue", zorder=2)
    for s, (e, n, _) in zip(array.stations, enu):
        ant = s.antenna.name if s.antenna else "-"
        ax.annotate(f"{s.name}\n{ant}", (e * f, n * f), textcoords="offset points",
                    xytext=(6, -14), fontsize=8)
    ax.set_xlabel(f"east  [{unit}]")
    ax.set_ylabel(f"north  [{unit}]")
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.25)
    ax.grid(alpha=0.3)
    ax.set_title(f"{array.name}: {len(array)} stations, site {array.site.name} "
                 f"({array.site.lat_deg:.5f}°, {array.site.lon_deg:.5f}°)", fontsize=9.5)
    return fig, ax
