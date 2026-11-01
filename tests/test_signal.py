# tests/test_signal.py -- tests of interferosim.observation.geometry and interferosim.signal
import math
from pathlib import Path

import numpy as np
import pytest

from interferosim.constants import C_LIGHT, HI_21CM
from interferosim.instrument import Antenna, Array, Receiver, Site, Station
from interferosim.observation import (delay_rate, direction, elevation_deg, geometric_delay,
                                      gha_deg, gmst_deg, sky_basis, uvw)
from interferosim.signal import (IQSimulator, StationErrors, correct_geometry, cross_correlate,
                                 fractional_delay, measure_visibilities, quantize)
from interferosim.sky import FieldCenter, PointSource, SkyModel

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs"


def _ew_array(b_m, antenna=None):
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    ant = antenna or Antenna.from_toml(CFG / "antennas" / "dish_3m.toml")
    return Array(site, [Station.from_enu("A", site, 0, 0, 0, ant),
                        Station.from_enu("B", site, b_m, 0, 0, ant)],
                 receiver=Receiver.from_toml(CFG / "receivers" / "b210.toml"))


# ---------------------------------------------------------------- geometry
def test_gmst_at_j2000():
    # 2000-01-01 12:00 UT: GMST = 280.46061837 deg
    assert float(gmst_deg(946728000.0)) == pytest.approx(280.46061837, abs=1e-6)


def test_sky_basis_orthonormal():
    s0, e_e, e_n = sky_basis(40.0, np.array([0.0, 37.0, 200.0]))
    for a in (s0, e_e, e_n):
        assert np.linalg.norm(a, axis=-1) == pytest.approx(1.0)
    assert np.einsum("ij,ij->i", s0, e_e) == pytest.approx(0.0, abs=1e-12)
    assert np.einsum("ij,ij->i", s0, e_n) == pytest.approx(0.0, abs=1e-12)


def test_delay_ew_baseline_transit_and_rate():
    arr = _ew_array(11_500.0)
    b = arr.baselines()[0]
    lon = arr.reference.lon_deg
    # at transit (LHA = 0) an east-west baseline has no delay, the rate is maximal
    g = float(gha_deg(0.0, lha_deg=0.0, lon_deg=lon))
    assert abs(float(geometric_delay(b, direction(30.0, g)))) < 1e-8
    rate = float(delay_rate(b, 30.0, g))
    assert abs(rate) == pytest.approx(7.292115e-5 * 11_500 * math.cos(math.radians(30)) / C_LIGHT, rel=2e-3)


def test_uvw_ew_baseline():
    arr = _ew_array(1000.0)
    lam = C_LIGHT / HI_21CM
    g = float(gha_deg(0.0, lha_deg=0.0, lon_deg=arr.reference.lon_deg))
    u, v, w = uvw(arr.baselines()[0], 20.0, g, lam)
    assert u == pytest.approx(1000.0 / lam, rel=1e-3)
    assert abs(v) < 2.0 and abs(w) < 2.0


def test_elevation_at_transit():
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    g = float(gha_deg(0.0, lha_deg=0.0, lon_deg=site.lon_deg))
    el = float(elevation_deg(direction(40.0, g), site.lat_deg, site.lon_deg))
    assert el == pytest.approx(90.0 - abs(site.lat_deg - 40.0), abs=0.25)   # geodetic vs geocentric


# ---------------------------------------------------------------- signal
def test_fractional_delay_shifts():
    x = np.zeros(64, complex); x[10] = 1.0
    y = fractional_delay(x, 3.0)
    assert np.argmax(np.abs(y)) == 13


def test_quantize_one_bit_van_vleck():
    rng = np.random.default_rng(0)
    n, rho = 400_000, 0.1
    c = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
    a = math.sqrt(rho) * c + math.sqrt(1 - rho) * (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
    b = math.sqrt(rho) * c + math.sqrt(1 - rho) * (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)
    r = lambda x, y: abs(np.mean(x * np.conj(y))) / math.sqrt(np.mean(abs(x) ** 2) * np.mean(abs(y) ** 2))
    assert r(quantize(a, 1), quantize(b, 1)) / r(a, b) == pytest.approx(2 / math.pi, abs=0.03)
    assert r(quantize(a, 8), quantize(b, 8)) / r(a, b) == pytest.approx(1.0, abs=0.01)


def test_power_and_noise_free_visibility_match_model():
    sky = SkyModel.from_json(CFG / "skies" / "cyg_a.json")
    sim = IQSimulator(sky, _ew_array(300.0), sample_rate_sps=2e4, seed=4, receiver_noise=False)
    bl = sim.simulate(10.0)
    assert np.mean(np.abs(bl["A"].data) ** 2) == pytest.approx(sim.station_sky_power_jy(0), rel=0.02)
    r = measure_visibilities(bl, sim, 10.0)["A-B"]
    exp = sim.expected_visibility(0, 1, *r["uvw"][0, :2])
    assert abs(r["vis"][0] - exp) < 5 * r["sigma"] + 0.01 * abs(exp)


def test_noise_level():
    sky = SkyModel([PointSource("p", 0, 0, 100.0)], center=FieldCenter("p", 120.0, 50.0))
    sim = IQSimulator(sky, _ew_array(50.0), sample_rate_sps=2e4, seed=8)
    bl = sim.simulate(5.0)
    assert np.mean(np.abs(bl["A"].data) ** 2) == pytest.approx(
        sim.receiver.sefd_jy(sim.array.stations[0].antenna) + 100.0, rel=0.02)
    r = measure_visibilities(bl, sim, 0.05)["A-B"]
    assert np.std(r["vis"].real) == pytest.approx(r["sigma"], rel=0.2)


def test_lo_offset_and_its_correction():
    sky = SkyModel([PointSource("p", 0, 0, 1000.0)], center=FieldCenter("p", 120.0, 50.0))
    df = 50.0
    sim = IQSimulator(sky, _ew_array(30.0), sample_rate_sps=2e4, seed=9, receiver_noise=False,
                      errors={"B": StationErrors(lo_offset_hz=df)})
    bl = sim.simulate(2.0)
    raw = measure_visibilities(bl, sim, 2.0)["A-B"]["vis"][0]
    fixed = measure_visibilities(bl, sim, 2.0, lo_offsets_hz={"B": df})["A-B"]["vis"][0]
    assert abs(raw) < 50.0
    assert abs(fixed) == pytest.approx(1000.0, rel=0.05)


def test_station_errors_from_receiver_limits():
    rec = Receiver.from_toml(CFG / "receivers" / "b210_gpsdo.toml")
    e = StationErrors.from_receiver(rec, np.random.default_rng(1))
    assert abs(e.lo_offset_hz) <= 1e-9 * HI_21CM + 1e-9
    assert abs(e.clock_offset_s) <= 50e-9


def test_sigmf_written(tmp_path):
    sky = SkyModel([PointSource("p", 0, 0, 10.0)], center=FieldCenter("p", 120.0, 50.0))
    sim = IQSimulator(sky, _ew_array(10.0), sample_rate_sps=1e4)
    blk = sim.simulate(0.5)["A"]
    meta = blk.save_sigmf(tmp_path / "rec")
    data = np.fromfile(meta.with_suffix(".sigmf-data"), dtype="<c8")
    assert data.size == blk.data.size and np.allclose(data, blk.data)


def test_sky_without_position_rejected():
    sun = SkyModel.from_json(CFG / "skies" / "sun_21cm.json")
    with pytest.raises(ValueError):
        IQSimulator(sun, _ew_array(10.0))


# ---------------------------------------------------------------- correlator (0.6.0)
from interferosim.signal import (Accumulator, FXCorrelator, apply_fringe_solution,  # noqa: E402
                                 cross_spectrum, delay_rate_search, fringe_stop,
                                 integer_align, stokes_from_products)


def test_integer_align_and_fringe_stop():
    a = np.arange(10.0); b = np.arange(10.0) + 3
    a2, b2 = integer_align(a, b, 2)
    assert a2.size == b2.size == 8 and b2[0] == 5
    x = np.exp(2j * np.pi * 5.0 * np.arange(100) / 100.0)
    assert np.allclose(fringe_stop(x, 5.0, 100.0), 1.0)


def test_delay_rate_search_finds_offsets():
    rng = np.random.default_rng(1)
    n, fs = 2 ** 15, 1e4
    s = rng.standard_normal(n + 64) + 1j * rng.standard_normal(n + 64)
    d, r = 3.3, 0.7                                   # samples, Hz
    a = s[32:32 + n]
    b = fractional_delay(s, d)[32:32 + n] * np.exp(2j * np.pi * r * np.arange(n) / fs)
    sol = delay_rate_search(a, b, fs, max_delay_samples=10, max_rate_hz=2, n_rate=81)
    # quadratic peak interpolation (AAS6): sub-sample bias up to ~0.3 samples for white noise
    assert sol.delay_samples == pytest.approx(-d, abs=0.35)
    assert sol.rate_hz == pytest.approx(r, abs=0.05)
    aligned = apply_fringe_solution(b, sol, fs)
    vis, coh, nw = cross_spectrum(a, aligned, 64)
    assert np.mean(coh) > 0.9


def test_accumulator_units_and_stokes():
    rng = np.random.default_rng(2)
    x = (rng.standard_normal(64 * 400) + 1j * rng.standard_normal(64 * 400)) * math.sqrt(5.0 / 2)
    acc = Accumulator.create(64)
    acc.add(x, x)
    assert np.mean(acc.visibility().real / acc.window_power) == pytest.approx(5.0, rel=0.03)
    st = stokes_from_products(np.ones(3), np.ones(3), np.zeros(3), np.zeros(3))
    assert np.allclose(st["I"], 2) and np.allclose(st["Q"], 0)


def test_fx_correlator_with_errors():
    sky = SkyModel.from_json(CFG / "skies" / "cyg_a.json")
    sim = IQSimulator(sky, _ew_array(100.0), sample_rate_sps=5e4, seed=4,
                      errors={"B": StationErrors(lo_offset_hz=0.6, clock_offset_s=60e-6)})
    bl = sim.simulate(6.0)
    fx = FXCorrelator(nfft=32, t_int_s=2.0)
    sol = fx.search(bl, sim, "B", segment_s=3.0, max_delay_samples=8, max_rate_hz=2, n_rate=81)
    assert sol.delay_samples == pytest.approx(-60e-6 * 5e4, abs=0.2)
    assert sol.rate_hz == pytest.approx(0.6, abs=0.1)
    assert abs(sol.delay_samples + 3.0) < 0.15          # refined
    vd = fx.correlate(bl, sim, solutions={"B": sol})
    assert vd.n_chan == 32 and vd.n_rows == 3
    avg = vd.average_channels()
    exp = sim.expected_visibility(0, 1, *sim.uvw(3.0)[0, :2])
    assert abs(np.mean(avg.vis[:, 0]) - exp) < 5 * avg.sigma[0, 0] / math.sqrt(3) + 0.2 * abs(exp)
    raw = fx.correlate(bl, sim).average_channels()
    assert abs(np.mean(raw.vis[:, 0])) < 0.5 * abs(exp)


# ---------------------------------------------------------------- XF correlator (0.9.1)
def test_lag_correlation_peak_and_xf_equals_fx_without_delay():
    from interferosim.signal.correlate import Accumulator, fractional_delay, lag_correlation, xf_cross_spectrum
    rng = np.random.default_rng(2)
    n = 2 ** 16
    s = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    b = fractional_delay(s, 4.0)                                    # b[n] = s[n-4]
    lags, R = lag_correlation(s, b, 32)
    assert lags[np.argmax(np.abs(R))] == -4 and abs(R).max() > 0.95
    # no delay: XF and FX (rectangular window) give the same mean cross spectrum
    a2 = s + 0.5 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    xf = xf_cross_spectrum(a2, s, 32)
    acc = Accumulator.create(32, "rect"); acc.add(a2, s)
    fx = acc.visibility() / acc.window_power
    assert np.abs(xf).mean() == pytest.approx(np.abs(fx).mean(), rel=0.03)
    # a delay appears as the phase slope 2 pi f d of the XF spectrum
    f = np.fft.fftshift(np.fft.fftfreq(32))
    S = xf_cross_spectrum(s, fractional_delay(s, 2.5), 32)
    inner = np.abs(f) < 0.45
    slope = np.polyfit(f[inner], np.unwrap(np.angle(S[inner])), 1)[0] / (2 * np.pi)
    assert slope == pytest.approx(2.5, abs=0.02)


# ---------------------------------------------------------------- baseband (0.9.2)
def test_quantization_efficiency_known_values_and_measured():
    from interferosim.signal import quantization_efficiency, quantize
    assert quantization_efficiency(1) == pytest.approx(2 / np.pi, rel=1e-6)
    assert quantization_efficiency(2) == pytest.approx(0.881, abs=0.001)
    assert quantization_efficiency(3) == pytest.approx(0.963, abs=0.002)
    assert quantization_efficiency(8) > 0.999
    rng = np.random.default_rng(4)
    n, rho = 2 ** 19, 0.1
    s = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    mk = lambda: np.sqrt(rho) * s + np.sqrt(1 - rho) * (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    a, b = mk(), mk()
    c = lambda x, y: np.real(np.mean(x * np.conj(y)) / np.sqrt(np.mean(abs(x) ** 2) * np.mean(abs(y) ** 2)))
    assert c(quantize(a, 2), quantize(b, 2)) / c(a, b) == pytest.approx(0.881, abs=0.02)


def test_iq_imbalance_image_and_van_vleck():
    from interferosim.signal import image_rejection_db, iq_imbalance, van_vleck_1bit
    t = np.arange(8192)
    x = np.exp(2j * np.pi * (1024 / 8192) * t)
    X = np.abs(np.fft.fft(iq_imbalance(x, 0.5, 3.0))) ** 2
    assert 10 * np.log10(X[1024] / X[-1024]) == pytest.approx(image_rejection_db(0.5, 3.0), abs=0.2)
    assert image_rejection_db(0.0, 0.0) == float("inf")
    assert van_vleck_1bit(2 / np.pi * np.arcsin(0.6)) == pytest.approx(0.6)
