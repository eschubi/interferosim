# tests/test_deconvolution.py -- tests of interferosim.deconvolution (0.8.0)
"""CLEAN, clean beam, Wiener filter and reference images on synthetic data
(cyclic FFT convolution, no observation needed)."""

import numpy as np
import pytest

from interferosim.deconvolution import (CleanBeam, Deconvolver, HogbomClean, WienerFilter, reference_image,
                                        rms_error)
from interferosim.imaging import Image
from interferosim.sky import SkyModel
from interferosim.sky.components import PointSource

N, PIX = 128, 0.25


def conv(a, kernel):
    return np.real(np.fft.ifft2(np.fft.fft2(a) * np.fft.fft2(np.fft.ifftshift(kernel))))


def ring_beam(seed=0):
    """Dirty beam of a sparse, ring-like uv sampling (with Hermitian symmetry)."""
    rng = np.random.default_rng(seed)
    W = np.zeros((N, N))
    for ang in rng.uniform(0, np.pi, 60):
        for r in (14, 22):
            iu, iv = int(round(r * np.cos(ang))), int(round(r * np.sin(ang) * 0.8))
            W[iv % N, iu % N] += 1; W[-iv % N, -iu % N] += 1
    b = np.fft.fftshift(np.real(np.fft.ifft2(W)))
    return Image(b / b.max(), PIX, ""), W


def point_model(sources):
    m = np.zeros((N, N))
    for j, i, s in sources:
        m[j, i] += s
    return m


SOURCES = [(64, 64, 10.0), (70, 52, 6.0), (50, 80, 3.0)]


def test_clean_beam_fit_recovers_gaussian():
    cb0 = CleanBeam(2.0, 1.2, 30.0)
    cb = CleanBeam.fit(cb0.image(N, PIX))
    assert cb.bmaj_arcsec == pytest.approx(2.0, rel=0.02)
    assert cb.bmin_arcsec == pytest.approx(1.2, rel=0.02)
    assert cb.pa_deg == pytest.approx(30.0, abs=1.0)


def test_clean_beam_convolve_keeps_point_flux_as_peak():
    cb = CleanBeam(1.5, 1.5)
    out = cb.convolve(point_model([(64, 64, 7.0)]), PIX)
    assert out.max() == pytest.approx(7.0, rel=1e-6)


def test_hogbom_recovers_point_sources():
    beam, _ = ring_beam()
    dirty = Image(conv(point_model(SOURCES), beam.data), PIX, "Jy/beam")
    r = HogbomClean(gain=0.2, niter=3000, threshold_jy=0.01).deconvolve(dirty, beam)
    a = dirty.axis_arcsec()
    for j, i, s in SOURCES:
        assert r.flux_within(a[i], a[j], 0.6) == pytest.approx(s, rel=0.03)
    assert r.residual.rms() < 0.01
    assert r.total_flux_jy == pytest.approx(19.0, rel=0.03)
    # history: residual peak decreases overall, flux grows to the total
    assert r.history["peak"][-1] < r.history["peak"][0]
    assert r.history["flux"][-1] == pytest.approx(r.total_flux_jy)


def test_hogbom_restored_shows_flux_at_positions():
    beam, _ = ring_beam()
    dirty = Image(conv(point_model(SOURCES), beam.data), PIX, "Jy/beam")
    r = HogbomClean(gain=0.1, niter=3000, threshold_jy=0.01).deconvolve(dirty, beam)
    a = dirty.axis_arcsec()
    j, i, s = SOURCES[0]
    assert r.restored.value_at(a[i], a[j]) == pytest.approx(s, rel=0.05)


def test_hogbom_window_and_positive():
    beam, _ = ring_beam()
    dirty = Image(conv(point_model(SOURCES), beam.data), PIX, "Jy/beam")
    a = dirty.axis_arcsec()
    r = HogbomClean(gain=0.2, niter=2000, window=[(a[64], a[64], 0.6)]).deconvolve(dirty, beam)
    for l, m, f in r.components():
        assert np.hypot(l - a[64], m - a[64]) <= 0.6 + 1e-9
    rp = HogbomClean(gain=0.2, niter=2000, positive=True, threshold_jy=0.01).deconvolve(dirty, beam)
    assert np.all(rp.model.data >= 0)


def test_hogbom_input_checks():
    beam, _ = ring_beam()
    with pytest.raises(ValueError):
        HogbomClean().deconvolve(Image(np.zeros((N // 2, N // 2)), PIX), beam)
    with pytest.raises(ValueError):
        HogbomClean(gain=0.0).deconvolve(Image(np.zeros((N, N)), PIX), beam)


def test_wiener_filled_aperture_sharpens():
    wide = CleanBeam(2.0, 2.0)
    narrow = CleanBeam(1.2, 1.2)
    model = point_model([(64, 64, 5.0), (64, 72, 5.0)])           # 2" apart, unresolved by 2"
    dirty = Image(wide.convolve(model, PIX), PIX, "Jy/beam")
    beam = wide.image(N, PIX)
    out = WienerFilter(1e-10).apply(dirty, beam, target=narrow)
    a = dirty.axis_arcsec()
    assert out.value_at(a[64], a[64]) == pytest.approx(5.0, rel=0.05)
    assert out.value_at(a[68], a[64]) < 0.4 * dirty.value_at(a[68], a[64])   # dip between the sources


def test_wiener_cannot_fill_unsampled_cells():
    beam, W = ring_beam()
    dirty = Image(conv(point_model(SOURCES), beam.data), PIX, "Jy/beam")
    out = WienerFilter(1e-3).apply(dirty, beam)
    F = np.fft.fft2(out.data)
    assert np.max(np.abs(F[W == 0])) < 1e-6 * np.max(np.abs(F))


def test_reference_image_and_rms_error():
    sky = SkyModel([PointSource("P", 0.0, 0.0, 4.0)])
    like = Image(np.zeros((N, N)), PIX)
    ref = reference_image(sky, like, clean_beam=CleanBeam(1.0, 1.0))
    assert ref.peak()[0] == pytest.approx(4.0, rel=1e-6)
    assert ref.peak()[1:] == (0.0, 0.0)
    assert rms_error(ref, ref) == 0.0
    assert rms_error(ref, like, radius_arcsec=0.1) == pytest.approx(4.0)


def test_deconvolver_interface_and_wiener_result():
    beam, _ = ring_beam()
    dirty = Image(conv(point_model(SOURCES), beam.data), PIX, "Jy/beam")
    for d in (HogbomClean(gain=0.2, niter=200), WienerFilter(1e-3)):
        assert isinstance(d, Deconvolver)
        r = d.run(dirty, beam)
        assert r.restored.data.shape == dirty.data.shape
        assert r.model.unit == "Jy/pixel" and r.restored.unit == "Jy/beam"
    # Wiener: model * beam reproduces the dirty image on the measured cells -> small residual
    rw = WienerFilter(1e-6).deconvolve(dirty, beam)
    assert rw.residual.rms() < 0.05 * dirty.rms()
