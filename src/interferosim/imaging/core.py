# src/interferosim/imaging/core.py -- gridding, weighting, dirty image and beam (interferosim.imaging.core)
"""Imaging of visibilities: from a :class:`~interferosim.visdata.VisData` set
to the dirty image and the dirty beam.

* :class:`Image` -- a square map with pixel size in arcsec (axes: offset
  from the image centre, l to the east, m to the north).
* :class:`Gridder` -- nearest-cell gridding on an ``npix`` x ``npix`` uv-grid
  with the cell 1 / (npix * cell); every visibility is entered together with
  its Hermitian conjugate V(-u, -v) = V*(u, v); all channels of the data set
  are gridded at their own u, v (multi-frequency synthesis) or the channel
  average.
* Weighting (strategy pattern, chapter "Basic Design"):
  :class:`NaturalWeighting` (data weights 1/sigma^2), :class:`UniformWeighting`
  (each occupied uv cell counts once), :class:`WienerWeighting` (W / (W^2 + K)
  per cell -- notebook A: for a sparsely sampled interferometer the Wiener
  filter only re-weights the measured cells).
* :class:`Imager` -- dirty image and dirty beam, normalised so that the beam
  peaks at 1 and a point source shows its flux density (Jy/beam); optional
  shift of the phase centre to the image centre; since 0.9.0
  :meth:`Imager.channel_images` (one dirty image per channel, a cube for
  spectral lines).

Origin: notebook A (``Fall B``: gridding with ``np.add.at``, ``image_from``,
``wiener_weights``); conventions as ``SkyModel.visibility``:
V(u, v) = sum S exp(-2 pi i (u l + v m))  =>  I(l, m) = sum V exp(+2 pi i (u l + v m)).

Version 0.7.0 (2026-10-31): first implementation.
Version 0.9.0 (2026-11-01): Imager.channel_images.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..visdata import VisData

__all__ = ["Image", "Gridder", "Weighting", "NaturalWeighting", "UniformWeighting",
           "WienerWeighting", "Imager"]

ARCSEC = math.pi / (180.0 * 3600.0)


# ======================================================================
@dataclass
class Image:
    """Square map, data[j, i] with i along l (east) and j along m (north)."""

    data: np.ndarray
    pix_arcsec: float
    unit: str = "Jy/beam"
    center_deg: Tuple[float, float] = (0.0, 0.0)     # offset of the image centre from the field centre
    name: str = ""
    meta: Dict[str, object] = field(default_factory=dict)

    @property
    def npix(self) -> int:
        return self.data.shape[0]

    def axis_arcsec(self) -> np.ndarray:
        """Pixel centres along l and m relative to the image centre [arcsec]."""
        return (np.arange(self.npix) - self.npix // 2) * self.pix_arcsec

    def extent_arcsec(self) -> Tuple[float, float, float, float]:
        a = self.axis_arcsec()
        h = 0.5 * self.pix_arcsec
        return (a[0] - h, a[-1] + h, a[0] - h, a[-1] + h)

    def value_at(self, l_arcsec: float, m_arcsec: float) -> float:
        """Pixel value nearest to the offset (l, m) from the image centre."""
        a = self.axis_arcsec()
        return float(self.data[int(np.argmin(np.abs(a - m_arcsec))), int(np.argmin(np.abs(a - l_arcsec)))])

    def peak(self) -> Tuple[float, float, float]:
        """(value, l, m) of the maximum."""
        j, i = np.unravel_index(int(np.argmax(self.data)), self.data.shape)
        a = self.axis_arcsec()
        return float(self.data[j, i]), float(a[i]), float(a[j])

    def rms(self, exclude_radius_arcsec: float = 0.0,
            around: Optional[Tuple[float, float]] = None) -> float:
        """rms of the map outside a circle (e.g. around the source)."""
        a = self.axis_arcsec()
        L, M = np.meshgrid(a, a)
        l0, m0 = around or (0.0, 0.0)
        mask = np.hypot(L - l0, M - m0) > exclude_radius_arcsec
        return float(np.sqrt(np.mean(self.data[mask] ** 2)))

    def beam_fwhm_arcsec(self) -> Tuple[float, float]:
        """FWHM of the main lobe along l and m (for a beam image: half-power
        crossings through the centre)."""
        c = self.npix // 2
        a = self.axis_arcsec()
        out = []
        for prof in (self.data[c, c:], self.data[c:, c]):
            k = int(np.argmax(prof < 0.5 * prof[0])) if np.any(prof < 0.5 * prof[0]) else len(prof) - 1
            if k == 0:
                out.append(float("nan")); continue
            x = np.interp(0.5 * prof[0], [prof[k], prof[k - 1]], [a[c + k], a[c + k - 1]])
            out.append(2.0 * float(x))
        return out[0], out[1]


# ======================================================================
class Weighting(ABC):
    """Interface: weights per uv cell from the gridded data weights W."""

    name: str = ""

    @abstractmethod
    def weights(self, W: np.ndarray) -> np.ndarray:
        """Factor per cell applied to the gridded (weighted) visibilities."""


class NaturalWeighting(Weighting):
    """Every visibility with its data weight 1/sigma^2 (best sensitivity,
    widest beam)."""

    name = "natural"

    def weights(self, W: np.ndarray) -> np.ndarray:
        return np.where(W > 0, 1.0, 0.0)


class UniformWeighting(Weighting):
    """Every occupied uv cell counts once (sharpest beam, more noise)."""

    name = "uniform"

    def weights(self, W: np.ndarray) -> np.ndarray:
        return np.where(W > 0, 1.0 / np.maximum(W, 1e-300), 0.0)


@dataclass
class WienerWeighting(Weighting):
    """Factor W / (W^2 + K) per cell (notebook A, ``wiener_weights``): K -> 0
    gives uniform weighting, K = <W>^2 a compromise between uniform and
    natural, larger K weights the cells increasingly with W (well-sampled cells
    dominate). ``k`` is given relative to the squared mean occupied-cell weight
    (k = 1: K = <W>^2)."""

    k: float = 1.0
    name: str = "wiener"

    def weights(self, W: np.ndarray) -> np.ndarray:
        occ = W > 0
        if not occ.any():
            return np.zeros_like(W)
        K = self.k * float(W[occ].mean()) ** 2
        return np.where(occ, 1.0 / (W + K / np.maximum(W, 1e-300)), 0.0)


# ======================================================================
@dataclass
class Gridder:
    """Nearest-cell gridder on an npix x npix grid with pixel ``cell_arcsec``."""

    npix: int = 512
    cell_arcsec: float = 0.5

    @property
    def du(self) -> float:
        """uv cell in wavelengths."""
        return 1.0 / (self.npix * self.cell_arcsec * ARCSEC)

    def grid(self, vd: VisData, center_deg: Tuple[float, float] = (0.0, 0.0),
             average_channels: bool = False) -> Tuple[np.ndarray, np.ndarray, int]:
        """Return (V, W, n_dropped): sum of w V and sum of w per cell (with the
        Hermitian conjugates), and the number of visibilities outside the grid.
        ``center_deg`` shifts the phase centre to the offset (l0, m0) of the
        field centre, so that this point lands in the image centre."""
        d = vd.average_channels() if average_channels else vd
        n = self.npix
        V = np.zeros((n, n), complex)
        W = np.zeros((n, n))
        l0, m0 = (math.radians(x) for x in center_deg)
        dropped = 0
        for k in range(d.n_chan):
            s = d.freq_hz[k] / d.freq_hz[0]
            u, v = d.u * s, d.v * s
            vis = d.vis[:, k] * np.exp(2j * np.pi * (u * l0 + v * m0))
            w = d.weight[:, k]
            ok = np.isfinite(vis) & (w > 0)
            iu = np.round(u / self.du).astype(int)
            iv = np.round(v / self.du).astype(int)
            inside = ok & (np.abs(iu) < n // 2) & (np.abs(iv) < n // 2)
            dropped += int(np.count_nonzero(ok & ~inside))
            iu, iv, vis, w = iu[inside], iv[inside], vis[inside], w[inside]
            for sgn, val in ((1, vis), (-1, np.conj(vis))):
                np.add.at(V, ((sgn * iv) % n, (sgn * iu) % n), w * val)
                np.add.at(W, ((sgn * iv) % n, (sgn * iu) % n), w)
        return V, W, dropped


# ======================================================================
@dataclass
class Imager:
    """Dirty image and dirty beam with a gridder and a weighting."""

    gridder: Gridder = field(default_factory=Gridder)
    weighting: Weighting = field(default_factory=NaturalWeighting)

    def _images(self, vd: VisData, center_deg=(0.0, 0.0), average_channels=False):
        V, W, dropped = self.gridder.grid(vd, center_deg, average_channels)
        f = self.weighting.weights(W)
        img = np.fft.fftshift(np.real(np.fft.ifft2(V * f)))
        beam = np.fft.fftshift(np.real(np.fft.ifft2(W * f)))
        p = float(beam.max()) or 1.0
        meta = {"weighting": self.weighting.name, "n_cells": int(np.count_nonzero(W)),
                "dropped": dropped, "du_lambda": self.gridder.du}
        return img / p, beam / p, meta

    def dirty_image(self, vd: VisData, center_deg: Tuple[float, float] = (0.0, 0.0),
                    average_channels: bool = False) -> Image:
        img, _, meta = self._images(vd, center_deg, average_channels)
        return Image(img, self.gridder.cell_arcsec, "Jy/beam", tuple(center_deg),
                     name=f"dirty image ({self.weighting.name})", meta=meta)

    def dirty_beam(self, vd: VisData, center_deg: Tuple[float, float] = (0.0, 0.0),
                   average_channels: bool = False) -> Image:
        _, beam, meta = self._images(vd, center_deg, average_channels)
        return Image(beam, self.gridder.cell_arcsec, "", tuple(center_deg),
                     name=f"dirty beam ({self.weighting.name})", meta=meta)

    def both(self, vd: VisData, center_deg: Tuple[float, float] = (0.0, 0.0),
             average_channels: bool = False) -> Tuple[Image, Image]:
        """(dirty image, dirty beam) in one gridding pass."""
        img, beam, meta = self._images(vd, center_deg, average_channels)
        return (Image(img, self.gridder.cell_arcsec, "Jy/beam", tuple(center_deg),
                      name=f"dirty image ({self.weighting.name})", meta=meta),
                Image(beam, self.gridder.cell_arcsec, "", tuple(center_deg),
                      name=f"dirty beam ({self.weighting.name})", meta=meta))

    def channel_images(self, vd: VisData, channels: Optional[Sequence[int]] = None,
                       center_deg: Tuple[float, float] = (0.0, 0.0)) -> List[Image]:
        """One dirty image per channel (all channels by default): the image
        cube of a spectral-line observation. ``meta["freq_hz"]`` holds the
        channel frequency."""
        out = []
        for k in (range(vd.n_chan) if channels is None else channels):
            img = self.dirty_image(vd.channel(int(k)), center_deg)
            img.meta = dict(img.meta, channel=int(k), freq_hz=float(vd.freq_hz[int(k)]))
            img.name = f"channel {int(k)}"
            out.append(img)
        return out
