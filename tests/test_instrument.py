# tests/test_instrument.py -- tests of interferosim.instrument and interferosim.plotting.instrument
import math
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from interferosim.constants import C_LIGHT, HI_21CM, OH_1665  # noqa: E402
from interferosim.instrument import (Antenna, Array, Horn, ParabolicDish,  # noqa: E402
                                     Receiver, Site, Station, Yagi, antenna_from_dict,
                                     ecef_to_geodetic, geodetic_to_ecef, otf_circular)

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs"


# ---------------------------------------------------------------- antennas
@pytest.mark.parametrize("name, cls", [("horn_15dbi", Horn), ("dish_3m", ParabolicDish),
                                       ("yagi", Yagi)])
def test_antenna_configs_load(name, cls):
    a = Antenna.from_toml(CFG / "antennas" / f"{name}.toml")
    assert isinstance(a, cls)
    assert all(math.isfinite(v) and v > 0 for v in a.summary().values())


def test_horn_matches_notebook_b():
    h = Antenna.from_toml(CFG / "antennas" / "horn_15dbi.toml")
    assert h.effective_area() == pytest.approx(0.112, rel=1e-2)       # aas6sim a_eff_m2
    assert h.hpbw_deg() == pytest.approx(
        math.degrees(math.sqrt(4 * math.pi / 10**1.5 / 1.133)), rel=1e-3)


def test_dish_matches_aas6sim():
    d = ParabolicDish("d", 3.0, 0.60, 1.20)
    lam = C_LIGHT / HI_21CM
    assert d.effective_area() == pytest.approx(0.60 * math.pi * 9 / 4)
    assert d.hpbw_deg() == pytest.approx(math.degrees(1.2 * lam / 3.0))
    assert d.u_max() == pytest.approx(3.0 / lam)
    assert d.u_cut_eff() == pytest.approx(3.0 / lam * 1.02 / 1.2)
    # gain and effective area consistent: A_eff = G lambda^2 / 4 pi
    assert d.gain() * lam**2 / (4 * math.pi) == pytest.approx(d.effective_area())


@pytest.mark.parametrize("ant", [Horn("h", 15.0), Yagi("y", 13, 0.85, 15.0),
                                 ParabolicDish("d", 3.0, 0.6, 1.2)])
def test_beam_half_power_at_hpbw(ant):
    h = ant.hpbw_deg()
    assert ant.beam(0.0) == pytest.approx(1.0)
    assert ant.beam(h / 2) == pytest.approx(0.5, abs=0.02)


def test_horn_fixed_aperture_scaling():
    h = Horn("h", 15.0, design_freq_hz=HI_21CM)
    assert h.gain_dbi(2 * HI_21CM) == pytest.approx(15.0 + 20 * math.log10(2))
    assert Horn("h", 15.0).gain_dbi(OH_1665) == 15.0


def test_otf_limits():
    assert otf_circular(0.0, 10.0) == pytest.approx(1.0)
    assert otf_circular(10.0, 10.0) == pytest.approx(0.0)
    assert otf_circular(20.0, 10.0) == 0.0


def test_unknown_antenna_type():
    with pytest.raises(ValueError):
        antenna_from_dict({"type": "helix"})


# ---------------------------------------------------------------- receiver
def test_b210_and_sefd():
    r = Receiver.from_toml(CFG / "receivers" / "b210.toml")
    assert r.rx_channels == 2 and r.t_sys_k == 120.0
    d = ParabolicDish("d", 3.0, 0.6, 1.2)
    sefd = 2 * 1.380649e-23 * 120.0 / d.effective_area() / 1e-26
    assert r.sefd_jy(d) == pytest.approx(sefd)
    assert r.baseline_sigma_jy(d, d, 1.0) == pytest.approx(sefd / math.sqrt(2 * 2e6))
    assert r.single_dish_sigma_jy(d, 1.0) == pytest.approx(sefd / math.sqrt(2e6))


# ---------------------------------------------------------------- geometry
def test_ecef_roundtrip():
    x = geodetic_to_ecef(48.088056, 15.756111, 470.0)
    lat, lon, h = ecef_to_geodetic(*x)
    assert (lat, lon) == pytest.approx((48.088056, 15.756111), abs=1e-9)
    assert h == pytest.approx(470.0, abs=1e-4)


def test_station_from_enu():
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    s = Station.from_enu("x", site, 100.0, -50.0, 3.0)
    assert s.enu_m(site) == pytest.approx([100.0, -50.0, 3.0], abs=1e-6)


def test_aas6_stations_array():
    arr = Array.from_toml(CFG / "sites" / "aas6_sim_stations.toml")
    assert arr.baseline_names() == ["A-B"]
    assert arr.site.name == "Wilhelmsburg"
    assert isinstance(arr.stations[0].antenna, ParabolicDish)
    assert arr.receiver.name == "Ettus USRP B210"
    e, n, u = arr.baselines_enu()[0]
    assert e == pytest.approx(11_093, abs=5)          # about 11.1 km east
    assert abs(n) < 20 and abs(u) < 20


def test_two_horns_array():
    arr = Array.from_toml(CFG / "sites" / "two_horns_10m.toml")
    assert arr.baseline_lengths_m()[0] == pytest.approx(10.0, abs=1e-6)
    assert arr.baselines_enu()[0] == pytest.approx([10.0, 0.0, 0.0], abs=1e-6)
    assert all(isinstance(s.antenna, Horn) for s in arr.stations)


def test_antenna_override_and_baseline_count():
    site = Site("s", 48.0, 15.0, 0.0)
    y = Yagi("y", 13, 0.85, 15.0)
    arr = Array(site, [Station.from_enu(f"S{k}", site, 10.0 * k, 0.0, 0.0, y) for k in range(4)])
    assert len(arr.baseline_pairs()) == 6
    assert np.sort(arr.baseline_lengths_m()) == pytest.approx([10, 10, 10, 20, 20, 30], abs=1e-6)


# ---------------------------------------------------------------- plotting
def test_instrument_plots_run():
    from interferosim.plotting import plot_array, plot_beams
    ants = [Antenna.from_toml(p) for p in sorted((CFG / "antennas").glob("*.toml"))]
    fig, axes = plot_beams(ants)
    assert len(axes) == 2
    fig, ax = plot_array(Array.from_toml(CFG / "sites" / "aas6_sim_stations.toml"))
    assert "km" in ax.get_xlabel()


# ---------------------------------------------------------------- B210 device profile (0.3.2)
@pytest.mark.parametrize("name", ["b210", "b210_gpsdo", "b210_aas6_sim"])
def test_b210_profiles_valid(name):
    r = Receiver.from_toml(CFG / "receivers" / f"{name}.toml")
    assert r.validate() == []
    assert r.adc_bits == 12 and r.freq_max_hz == 6.0e9      # inherited via base
    assert r.extra == {}


def test_b210_profile_inheritance():
    base = Receiver.from_toml(CFG / "receivers" / "b210.toml")
    gps = Receiver.from_toml(CFG / "receivers" / "b210_gpsdo.toml")
    sim = Receiver.from_toml(CFG / "receivers" / "b210_aas6_sim.toml")
    assert gps.clock_source == "gpsdo" and gps.rx_channels == 1
    assert gps.sample_rate_sps == base.sample_rate_sps
    assert sim.center_freq_hz == pytest.approx(1665.402e6)
    assert sim.sample_rate_sps == 2.0e5 and sim.rx_gain_db == 40.0


def test_b210_derived_quantities():
    r = Receiver.from_toml(CFG / "receivers" / "b210.toml")
    assert r.data_rate_bytes_per_s() == pytest.approx(2.0e6 * 4 * 2)     # sc16, two channels
    assert r.receiver_noise_temperature_k() == pytest.approx(290 * (10**0.8 - 1))
    assert r.lo_offset_hz() == pytest.approx(2e-6 * HI_21CM)
    assert r.lo_offset_hz(clock_source="gpsdo") == pytest.approx(1e-9 * HI_21CM)
    t1 = r.phase_coherence_time_s(clock_source="gpsdo")
    assert t1 == pytest.approx(1 / (2 * math.pi * 1e-9 * HI_21CM))


def test_receiver_validate_flags_limits(tmp_path):
    r = Receiver(sample_rate_sps=40e6, bandwidth_hz=50e6, rx_channels=2,
                 max_sample_rate_2ch_sps=30.72e6, max_bandwidth_2ch_hz=30.72e6,
                 center_freq_hz=10e9, freq_min_hz=70e6, freq_max_hz=6e9,
                 rx_gain_db=90.0, rx_gain_range_db=[0.0, 76.0], noise_figure_db=8.0)
    w = r.validate()
    assert len(w) == 6
    p = tmp_path / "a.toml"
    p.write_text('base = "a"\n', encoding="utf-8")
    with pytest.raises(ValueError):
        Receiver.from_toml(p)
