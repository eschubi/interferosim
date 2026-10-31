# tests/test_spectral.py -- spectra, Doppler correction, multichannel visibilities (0.9.0)
"""Spectra of sky components, LSR correction, SpectralVisibilitySimulator,
regridding to the LSR and channel images."""

from pathlib import Path

import numpy as np
import pytest

from interferosim.constants import OH_1665
from interferosim.imaging import Gridder, Imager
from interferosim.instrument import Antenna, Array, Receiver, Site, Station
from interferosim.observation import (Observation, SpectralVisibilitySimulator, observer_velocity_kms,
                                      regrid_to_lsr, topocentric_frequency_hz, v_lsr_correction_kms)
from interferosim.sky import (AbsorptionSpectrum, FieldCenter, GaussianLineSpectrum, PointSource,
                              PowerLawSpectrum, SkyModel, frequency_from_velocity,
                              radio_velocity_kms, spectrum_from_dict)

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs"


# ---------------------------------------------------------------- spectra
def test_velocity_conversion_roundtrip():
    f = frequency_from_velocity(-44.0, OH_1665)
    assert radio_velocity_kms(f, OH_1665) == pytest.approx(-44.0)
    assert f > OH_1665                                    # approaching -> higher frequency


def test_spectrum_shapes_and_json_roundtrip():
    line = GaussianLineSpectrum(OH_1665, -44.0, 0.4)
    f0 = frequency_from_velocity(-44.0, OH_1665)
    assert line.shape(f0) == pytest.approx(1.0)
    assert line.shape(frequency_from_velocity(-43.8, OH_1665)) == pytest.approx(0.5)   # half width
    pl = PowerLawSpectrum(-0.7, 1.4e9)
    assert pl.shape(2.8e9) == pytest.approx(2 ** -0.7)
    ab = AbsorptionSpectrum(1.42040575e9, [(0.0, 5.0, 0.5)])
    assert ab.shape(1.42040575e9) == pytest.approx(np.exp(-0.5))
    for sp in (line, pl, ab):
        assert spectrum_from_dict(sp.to_dict()) == sp


def test_skymodel_spectra_from_json_and_at_frequency():
    sky = SkyModel.from_json(CFG / "skies" / "w3oh_masers_lines.json")
    assert set(sky.spectra) == {"Maser M1", "Maser M2", "Maser M3"}
    f1 = frequency_from_velocity(sky.spectra["Maser M1"].v_kms, OH_1665)
    assert float(sky.flux_at("Maser M1", f1)) == pytest.approx(70.0)
    snap = sky.at_frequency(f1)
    assert snap["Maser M1"].flux_jy == pytest.approx(70.0) and snap["Maser M3"].flux_jy < 1e-6
    assert not snap.spectra
    back = SkyModel.from_dict(sky.to_dict())
    assert back.spectra == sky.spectra
    # components without a spectrum stay flat
    plain = SkyModel.from_json(CFG / "skies" / "w3oh_masers.json")
    assert float(plain.flux_at("Maser M1", 1.0e9)) == pytest.approx(70.0)


def test_visibility_spectrum_scales_uv_with_frequency():
    sky = SkyModel([PointSource("P", 1e-3, 0.0, 2.0)])
    sky.set_spectrum("P", PowerLawSpectrum(-1.0, 1.0e9))
    V = sky.visibility_spectrum([1000.0], [0.0], [1.0e9, 2.0e9], ref_hz=1.0e9)
    assert abs(V[0, 1]) == pytest.approx(1.0)
    assert abs(V[0, 0]) == pytest.approx(2.0)
    # the phase -2 pi u l doubles with the frequency (u in wavelengths scales with f)
    assert np.angle(V[0, 1]) / np.angle(V[0, 0]) == pytest.approx(2.0, rel=1e-9)


# ---------------------------------------------------------------- Doppler
def test_observer_velocity_magnitudes():
    t = 61345.5 + np.linspace(0, 365, 50)
    orb = observer_velocity_kms(t, 48.11, 15.76, 280, components=("orbit",))
    assert np.linalg.norm(orb, axis=1).min() == pytest.approx(29.3, abs=0.3)
    assert np.linalg.norm(orb, axis=1).max() == pytest.approx(30.3, abs=0.3)
    rot = observer_velocity_kms(t, 48.11, 15.76, 280, components=("rotation",))
    assert np.linalg.norm(rot, axis=1) == pytest.approx(0.465 * np.cos(np.radians(48.11)), rel=0.01)
    sun = observer_velocity_kms(61345.5, 0, 0, 0, components=("sun",))
    assert np.linalg.norm(sun) == pytest.approx(20.0)


def test_lsr_correction_annual_and_inverse():
    ra, dec = 299.868, 40.733                               # Cyg A
    t = 61345.5 + np.arange(0, 366, 5.0)
    vc = v_lsr_correction_kms(ra, dec, t, 48.11, 15.76, 280)
    # orbital amplitude projected on a source at ecliptic latitude ~ +59 deg: ~15 km/s
    assert 12 < 0.5 * (vc.max() - vc.min()) < 18
    f = topocentric_frequency_hz(OH_1665, -44.0, ra, dec, 61345.5, 48.11, 15.76, 280)
    from interferosim.observation import lsr_velocity_from_topocentric
    assert float(lsr_velocity_from_topocentric(f, OH_1665, ra, dec, 61345.5, 48.11, 15.76, 280)) \
        == pytest.approx(-44.0, abs=1e-9)


# ---------------------------------------------------------------- multichannel visibilities
def _array():
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    dish = Antenna.from_toml(CFG / "antennas" / "dish_3m.toml")
    rec = Receiver.from_toml(CFG / "receivers" / "b210.toml")
    return Array(site, [Station.from_enu("W", site, 0, 0, 0, dish),
                        Station.from_enu("X", site, 11_500, 0, 0, dish)], receiver=rec)


def _maser_at_centre(v=-44.0):
    sky = SkyModel([PointSource("M", 0.0, 0.0, 50.0)], center=FieldCenter("W3(OH)", 36.7661, 61.8742))
    sky.set_spectrum("M", GaussianLineSpectrum(OH_1665, v, 0.5))
    return sky


def test_spectral_simulator_puts_line_at_lsr_velocity():
    arr, sky = _array(), _maser_at_centre()
    obs = Observation.from_sky(sky, arr, ha_range_h=(-4, 4), t_int_s=300, freq_hz=OH_1665, date="2026-11-01")
    sim = SpectralVisibilitySimulator(n_chan=128, bandwidth_hz=125e3, rest_hz=OH_1665, v_center_kms=-44.0,
                                      add_noise=False)
    vd = sim.run(sky, obs)
    assert vd.vis.shape == (vd.n_rows, 128)
    v = np.linspace(-46, -42, 81)
    rg = regrid_to_lsr(vd, arr, OH_1665, v)
    amp = np.abs(rg.vis)
    k = np.argmax(amp.mean(axis=0))
    assert v[k] == pytest.approx(-44.0, abs=0.06)
    # line stays put after regridding; linear interpolation of a 0.5 km/s line sampled every
    # 0.18 km/s loses at most ~8 % of the peak
    assert amp[:, k].min() > 0.9 * 50.0
    # without regridding the line drifts with the Earth's rotation (topocentric channels)
    peak_chan = np.argmax(np.abs(vd.vis), axis=1)
    assert peak_chan.max() - peak_chan.min() >= 1


def test_spectral_simulator_noise_per_channel():
    arr, sky = _array(), _maser_at_centre()
    obs = Observation.from_sky(sky, arr, ha_range_h=(-1, 1), t_int_s=60, freq_hz=OH_1665, date="2026-11-01")
    vd = SpectralVisibilitySimulator(n_chan=64, bandwidth_hz=64e3, rest_hz=OH_1665, v_center_kms=0.0,
                                     seed=3).run(sky, obs)
    sefd = arr.receiver.sefd_jy(arr.stations[0].antenna, OH_1665)
    expect = sefd / np.sqrt(2 * 1e3 * 60)
    assert vd.sigma[0, 0] == pytest.approx(expect, rel=1e-3)
    assert np.std(vd.vis.real) == pytest.approx(expect, rel=0.05)      # line far outside the band


def test_channel_images_find_line_channel():
    arr, sky = _array(), _maser_at_centre()
    obs = Observation.from_sky(sky, arr, ha_range_h=(-6, 6), t_int_s=300, freq_hz=OH_1665, date="2026-11-01")
    vd = SpectralVisibilitySimulator(n_chan=32, bandwidth_hz=62.5e3, rest_hz=OH_1665, v_center_kms=-44.0,
                                     add_noise=False).run(sky, obs)
    imgs = Imager(Gridder(128, 0.5)).channel_images(vd)
    peaks = [im.value_at(0, 0) for im in imgs]
    k = int(np.argmax(peaks))
    assert abs(k - 16) <= 1 and peaks[k] > 30.0 and imgs[k].meta["channel"] == k
