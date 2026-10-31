# src/interferosim/deconvolution/core.py -- CLEAN, clean beam, Wiener filter (interferosim.deconvolution.core)
"""Deconvolution of a dirty image with its dirty beam.

* :class:`Deconvolver` -- interface (basic design): ``deconvolve(dirty,
  beam) -> DeconvResult``, alias ``run``; realised by :class:`HogbomClean`
  and :class:`WienerFilter`.
* :class:`CleanBeam` -- elliptical Gaussian fitted to the main lobe of the
  dirty beam (major and minor axis FWHM, position angle east of north);
  kernel image and convolution of a model.
* :class:`HogbomClean` -- Hoegbom CLEAN (1974): find the peak of the
  residual, subtract the dirty beam shifted to it and scaled by the loop
  gain, store the component; stop at a threshold or after ``niter``
  iterations; optional CLEAN window and positivity.
* :class:`DeconvResult` -- component model (Jy/pixel), residual and restored
  image (model * clean beam + residual, Jy/beam), the clean beam and the
  convergence history.
* :class:`WienerFilter` -- linear deconvolution in the Fourier domain,
  ``W = S* / (|S|^2 + NSR)``, optionally to the resolution of a target beam;
  for a filled aperture (single dish) it restores resolution up to D/lambda,
  for a sparsely sampled interferometer it cannot fill the empty uv cells
  (notebook A). :meth:`WienerFilter.apply` returns the filtered image,
  :meth:`WienerFilter.deconvolve` a :class:`DeconvResult`.
* :func:`reference_image` -- sky model on the pixel grid of an image,
  convolved with a clean beam (the "truth" for comparisons), and
  :func:`rms_error` within a circle.

Origin: notebook A of the teaching material (``hogbom_fft``, clean beam from
the main lobe, ``conv_cbeam``, ``rms_err``). The dirty beam of
:class:`~interferosim.imaging.Imager` is periodic in the FFT grid, so
shifting it is a cyclic roll (implemented as a slice of the tiled beam).

Version 0.8.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..imaging import Image

__all__ = ["Deconvolver", "CleanBeam", "HogbomClean", "DeconvResult", "WienerFilter", "reference_image",
           "rms_error"]

FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))


# ======================================================================
def _grid_lm(npix: int, pix_arcsec: float) -> Tuple[np.ndarray, np.ndarray]:
    a = (np.arange(npix) - npix // 2) * pix_arcsec
    return np.meshgrid(a, a)


def _fft_convolve(data: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Cyclic convolution with a kernel centred at (npix//2, npix//2)."""
    return np.real(np.fft.ifft2(np.fft.fft2(data) * np.fft.fft2(np.fft.ifftshift(kernel))))


@dataclass
class CleanBeam:
    """Elliptical Gaussian, peak 1: FWHM of major and minor axis [arcsec],
    position angle of the major axis east of north [deg]."""

    bmaj_arcsec: float
    bmin_arcsec: float
    pa_deg: float = 0.0

    @classmethod
    def fit(cls, beam: Image, level: float = 0.35) -> "CleanBeam":
        """Least-squares fit of ln B = -(a l^2 + b l m + c m^2) to the pixels of
        the main lobe above ``level`` (connected region around the centre,
        bounded by twice the half-power width)."""
        n = beam.npix
        c0 = n // 2
        fx, fy = beam.beam_fwhm_arcsec()
        r = 1.2 * max(v for v in (fx, fy, 2 * beam.pix_arcsec) if np.isfinite(v))
        L, M = _grid_lm(n, beam.pix_arcsec)
        d = beam.data / beam.data[c0, c0]
        sel = (np.hypot(L, M) <= r) & (d > level)
        x, y, z = L[sel], M[sel], np.log(d[sel])
        A = np.column_stack([x * x, x * y, y * y])
        (a, b, c), *_ = np.linalg.lstsq(A, -z, rcond=None)
        Q = 2.0 * np.array([[a, b / 2.0], [b / 2.0, c]])          # inverse covariance
        lam, vec = np.linalg.eigh(Q)
        lam = np.maximum(lam, 1e-12)
        s_maj, s_min = 1.0 / math.sqrt(lam[0]), 1.0 / math.sqrt(lam[1])
        vx, vy = vec[:, 0]                                          # major axis direction (l, m)
        pa = math.degrees(math.atan2(vx, vy)) % 180.0               # east of north
        return cls(s_maj * FWHM_PER_SIGMA, s_min * FWHM_PER_SIGMA, pa)

    @property
    def area_pix(self) -> float:
        """Not used for the flux scale (peak-1 kernel); beam solid angle in arcsec^2."""
        return math.pi * self.bmaj_arcsec * self.bmin_arcsec / (4.0 * math.log(2.0))

    def kernel(self, npix: int, pix_arcsec: float) -> np.ndarray:
        """Gaussian of peak 1 centred at (npix//2, npix//2)."""
        L, M = _grid_lm(npix, pix_arcsec)
        pa = math.radians(self.pa_deg)
        # coordinates along the major (north-through-east) and minor axis
        x_maj = L * math.sin(pa) + M * math.cos(pa)
        x_min = L * math.cos(pa) - M * math.sin(pa)
        s1, s2 = self.bmaj_arcsec / FWHM_PER_SIGMA, self.bmin_arcsec / FWHM_PER_SIGMA
        return np.exp(-0.5 * ((x_maj / s1) ** 2 + (x_min / s2) ** 2))

    def image(self, npix: int, pix_arcsec: float) -> Image:
        return Image(self.kernel(npix, pix_arcsec), pix_arcsec, "", name="clean beam")

    def convolve(self, data: np.ndarray, pix_arcsec: float) -> np.ndarray:
        """Model in Jy/pixel -> image in Jy/beam (a point keeps its flux as peak)."""
        return _fft_convolve(data, self.kernel(data.shape[0], pix_arcsec))

    def __str__(self) -> str:
        return f'{self.bmaj_arcsec:.2f}" x {self.bmin_arcsec:.2f}", PA {self.pa_deg:.0f} deg'


# ======================================================================
@dataclass
class DeconvResult:
    """Result of a deconvolution."""

    model: Image                     # CLEAN components, Jy/pixel
    residual: Image                  # Jy/beam
    restored: Image                  # model * clean beam + residual, Jy/beam
    clean_beam: CleanBeam
    n_iter: int
    history: Dict[str, np.ndarray] = field(default_factory=dict)   # peak, flux per iteration
    meta: Dict[str, object] = field(default_factory=dict)

    @property
    def total_flux_jy(self) -> float:
        return float(self.model.data.sum())

    def components(self, min_flux_jy: float = 0.0) -> List[Tuple[float, float, float]]:
        """(l, m, flux) of all model pixels with ``|flux| > min_flux_jy``, brightest first."""
        a = self.model.axis_arcsec()
        j, i = np.nonzero(np.abs(self.model.data) > min_flux_jy)
        out = [(float(a[ii]), float(a[jj]), float(self.model.data[jj, ii])) for jj, ii in zip(j, i)]
        return sorted(out, key=lambda t: -abs(t[2]))

    def flux_within(self, l_arcsec: float, m_arcsec: float, radius_arcsec: float) -> float:
        """Sum of the model inside a circle (flux of a source broken into components)."""
        L, M = _grid_lm(self.model.npix, self.model.pix_arcsec)
        return float(self.model.data[np.hypot(L - l_arcsec, M - m_arcsec) <= radius_arcsec].sum())


# ======================================================================
class Deconvolver(ABC):
    """Interface of the deconvolution methods (strategy pattern, chapter
    "Basic Design")."""

    @abstractmethod
    def deconvolve(self, dirty: Image, beam: Image,
                   clean_beam: Optional[CleanBeam] = None) -> DeconvResult:
        """Deconvolve ``dirty`` with the dirty ``beam`` (same grid)."""

    def run(self, dirty: Image, beam: Image, clean_beam: Optional[CleanBeam] = None) -> DeconvResult:
        """Alias of :meth:`deconvolve` (name in the UML diagram)."""
        return self.deconvolve(dirty, beam, clean_beam)


def _result(dirty: Image, model: np.ndarray, residual: np.ndarray, cb: CleanBeam, n_iter: int,
            history: Dict[str, np.ndarray], meta: Dict[str, object], method: str) -> DeconvResult:
    restored = cb.convolve(model, dirty.pix_arcsec) + residual
    mk = lambda d, unit, name: Image(d, dirty.pix_arcsec, unit, dirty.center_deg, name=name, meta=meta)
    return DeconvResult(model=mk(model, "Jy/pixel", f"{method} model"),
                        residual=mk(residual, "Jy/beam", "residual"),
                        restored=mk(restored, "Jy/beam", f"restored image ({method})"),
                        clean_beam=cb, n_iter=n_iter, history=history, meta=meta)


# ======================================================================
WindowSpec = Union[None, np.ndarray, Sequence[Tuple[float, float, float]]]


def _window_mask(img: Image, window: WindowSpec) -> np.ndarray:
    if window is None:
        return np.ones(img.data.shape, bool)
    if isinstance(window, np.ndarray):
        if window.shape != img.data.shape:
            raise ValueError("window mask must have the shape of the image")
        return window.astype(bool)
    L, M = _grid_lm(img.npix, img.pix_arcsec)
    mask = np.zeros(img.data.shape, bool)
    for l0, m0, r in window:
        mask |= np.hypot(L - l0, M - m0) <= r
    return mask


@dataclass
class HogbomClean(Deconvolver):
    """Hoegbom CLEAN with loop gain ``gain``.

    Stops after ``niter`` iterations or when the residual peak inside the
    window falls below max(``threshold_jy``, ``frac`` x initial peak).
    ``window``: None (whole image), a boolean mask, or circles (l, m, r) in
    arcsec. ``positive``: search only positive peaks (masers, continuum
    without absorption)."""

    gain: float = 0.1
    niter: int = 1000
    threshold_jy: float = 0.0
    frac: float = 0.0
    positive: bool = False
    window: WindowSpec = None

    def deconvolve(self, dirty: Image, beam: Image, clean_beam: Optional[CleanBeam] = None) -> DeconvResult:
        if dirty.data.shape != beam.data.shape:
            raise ValueError("dirty image and dirty beam must have the same shape")
        if not 0.0 < self.gain <= 1.0:
            raise ValueError("gain must be in (0, 1]")
        n = dirty.npix
        c = n // 2
        bm = beam.data / beam.data[c, c]
        tiled = np.tile(bm, (2, 2))
        res = dirty.data.astype(float).copy()
        cc = np.zeros_like(res)
        mask = _window_mask(dirty, self.window)
        search = np.where(mask, res, 0.0)
        p0 = float(np.max(search) if self.positive else np.max(np.abs(search)))
        thr = max(self.threshold_jy, self.frac * p0)
        peaks, fluxes = [], []
        k = 0
        for k in range(1, self.niter + 1):
            s = np.where(mask, res, -np.inf if self.positive else 0.0)
            idx = int(np.argmax(s) if self.positive else np.argmax(np.abs(s)))
            j, i = divmod(idx, n)
            pk = res[j, i]
            if abs(pk) <= thr or (self.positive and pk <= 0):
                k -= 1
                break
            f = self.gain * pk
            cc[j, i] += f
            sj, si = (c - j) % n, (c - i) % n                     # start of the rolled beam in the tile
            res -= f * tiled[sj:sj + n, si:si + n]
            peaks.append(abs(pk)); fluxes.append(f)
        cb = clean_beam or CleanBeam.fit(beam)
        meta = dict(dirty.meta, gain=self.gain, threshold_jy=thr, initial_peak=p0)
        return _result(dirty, cc, res, cb, k, {"peak": np.array(peaks), "flux": np.cumsum(fluxes)},
                       meta, "CLEAN")


# ======================================================================
@dataclass
class WienerFilter(Deconvolver):
    """Wiener deconvolution ``W = S* / (|S|^2 + NSR)`` with the transfer function
    S = FFT(dirty beam). ``nsr`` is the noise-to-signal power ratio relative
    to ``max|S|^2`` (a number or an array on the FFT grid). With a ``target``
    clean beam the result has its resolution (Jy/beam); without, the result
    is the deconvolved map in Jy/pixel."""

    nsr: Union[float, np.ndarray] = 1e-2

    def transfer(self, beam: Image) -> np.ndarray:
        c = beam.npix // 2
        return np.fft.fft2(np.fft.ifftshift(beam.data / beam.data[c, c]))

    def apply(self, dirty: Image, beam: Image, target: Optional[CleanBeam] = None) -> Image:
        if dirty.data.shape != beam.data.shape:
            raise ValueError("dirty image and dirty beam must have the same shape")
        S = self.transfer(beam)
        p = float(np.max(np.abs(S)) ** 2)
        W = np.conj(S) / (np.abs(S) ** 2 + np.asarray(self.nsr) * p)
        T = 1.0 if target is None else np.fft.fft2(np.fft.ifftshift(target.kernel(dirty.npix, dirty.pix_arcsec)))
        out = np.real(np.fft.ifft2(np.fft.fft2(dirty.data) * W * T))
        unit = "Jy/pixel" if target is None else "Jy/beam"
        return Image(out, dirty.pix_arcsec, unit, dirty.center_deg,
                     name="Wiener filter", meta=dict(dirty.meta, nsr=self.nsr))

    def deconvolve(self, dirty: Image, beam: Image,
                   clean_beam: Optional[CleanBeam] = None) -> DeconvResult:
        """Deconvolved map as model (Jy/pixel), residual = dirty - model * beam,
        restored = model * clean beam + residual."""
        model = self.apply(dirty, beam).data
        c = beam.npix // 2
        res = dirty.data - _fft_convolve(model, beam.data / beam.data[c, c])
        cb = clean_beam or CleanBeam.fit(beam)
        return _result(dirty, model, res, cb, 1, {}, dict(dirty.meta, nsr=self.nsr), "Wiener")


# ======================================================================
def reference_image(sky, like: Image, center_deg: Tuple[float, float] = (0.0, 0.0),
                    clean_beam: Optional[CleanBeam] = None, antenna=None,
                    freq_hz: Optional[float] = None) -> Image:
    """Sky model rendered on the grid of ``like`` (image centre at the offset
    ``center_deg`` of the field centre), convolved with ``clean_beam``
    (Jy/beam) or left in Jy/pixel. With ``antenna`` and ``freq_hz`` every
    component is attenuated by the primary beam at its distance from the
    field centre (the apparent sky that the visibilities measure)."""
    a = like.axis_arcsec() / 3600.0
    x, y = a + center_deg[0], a + center_deg[1]
    pix = like.pix_arcsec / 3600.0
    if antenna is None:
        data = np.asarray(sky.render(x, y, pix), dtype=float)
    else:
        data = np.zeros((a.size, a.size))
        for comp in sky:
            part = np.zeros_like(data)
            comp.render(x, y, pix, part)
            theta = math.hypot(comp.dx_deg, comp.dy_deg)
            data += part * float(antenna.beam(np.asarray(theta), freq_hz))
    if clean_beam is not None:
        data = clean_beam.convolve(data, like.pix_arcsec)
    return Image(data, like.pix_arcsec, "Jy/beam" if clean_beam else "Jy/pixel", tuple(center_deg),
                 name="sky model" + (" * clean beam" if clean_beam else ""))


def rms_error(img: Image, ref: Image, radius_arcsec: Optional[float] = None,
              around: Tuple[float, float] = (0.0, 0.0)) -> float:
    """rms of img - ref inside a circle (whole image without radius)."""
    d = img.data - ref.data
    if radius_arcsec is None:
        return float(np.sqrt(np.mean(d ** 2)))
    L, M = _grid_lm(img.npix, img.pix_arcsec)
    m = np.hypot(L - around[0], M - around[1]) <= radius_arcsec
    return float(np.sqrt(np.mean(d[m] ** 2)))
