# tests/test_observation.py -- tests of interferosim.observation (obs, sun) and interferosim.visdata
import math
from pathlib import Path

import numpy as np
import pytest

from interferosim.constants import C_LIGHT
from interferosim.instrument import Array, Site, Station
from interferosim.observation import (Observation, TransitSimulator, UVCoverage,
                                      VisibilitySimulator, gmst_deg, mjd_from_iso, sun_radec_deg)
from interferosim.signal import IQSimulator, measure_visibilities
from interferosim.sky import FieldCenter, PointSource, SkyModel
from interferosim.visdata import VisData

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs"


def _aas6():
    return Array.from_toml(CFG / "sites" / "aas6_sim_stations.toml")


# ---------------------------------------------------------------- sun
def test_sun_position_equinox_and_solstice():
    ra, dec = sun_radec_deg(mjd_from_iso("2026-03-20T14:46:00"))      # March equinox 2026
    assert abs(float(dec)) < 0.05 and min(float(ra), 360 - float(ra)) < 0.1
    ra, dec = sun_radec_deg(mjd_from_iso("2026-06-21T08:24:00"))      # June solstice 2026
    assert float(dec) == pytest.approx(23.44, abs=0.05) and float(ra) == pytest.approx(90.0, abs=0.1)


# ---------------------------------------------------------------- observation
def test_transit_time_is_transit():
    obs = Observation(_aas6(), ra_deg=299.868, dec_deg=40.73, date="2026-11-01")
    t = obs.transit_mjd()
    lst = float(gmst_deg((t - 40587.0) * 86400.0)) + obs.array.reference.lon_deg
    assert (lst - 299.868 + 180) % 360 - 180 == pytest.approx(0.0, abs=1e-3)
    assert mjd_from_iso("2026-11-01") <= t < mjd_from_iso("2026-11-02")


def test_hour_angles_and_elevation_limit():
    obs = Observation(_aas6(), ra_deg=0.0, dec_deg=-30.0, ha_range_h=(-6, 6), t_int_s=600, min_elev_deg=5)
    ha = obs.hour_angles()
    assert ha[0] > -90 and ha[-1] < 90 and np.all(np.diff(ha) > 0)
    el = obs.elevation_deg()
    assert el.shape == (2, ha.size)
    # Dec -30 at latitude 48: maximum elevation about 12 deg
    assert el.max() == pytest.approx(90 - 48.088 - 30, abs=0.3)
    assert obs.visible().sum() < ha.size


def test_ew_uv_track_is_ellipse():
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    arr = Array(site, [Station.from_enu("A", site, 0, 0, 0), Station.from_enu("B", site, 1000, 0, 0)])
    dec = 50.0
    obs = Observation(arr, ra_deg=0.0, dec_deg=dec, ha_range_h=(-12, 12), t_int_s=600, min_elev_deg=-90)
    cov = UVCoverage.from_observation(obs)
    lam = C_LIGHT / obs.freq_hz
    u, v = cov.u[0], cov.v[0]
    b = 1000.0 / lam
    # east-west baseline: u^2 + (v / sin dec)^2 = b^2 (small corrections from the curved site frame)
    assert np.allclose(u ** 2 + (v / math.sin(math.radians(dec))) ** 2, b ** 2, rtol=2e-3)


def test_visibility_simulator_noise_free_matches_sky():
    sky = SkyModel.from_json(CFG / "skies" / "cyg_a.json")
    obs = Observation.from_sky(sky, _aas6(), ha_range_h=(-1, 1), t_int_s=600)
    vd = VisibilitySimulator(add_noise=False).run(sky, obs)
    exp = sky.visibility(vd.u, vd.v)       # beams ~1 for this small source
    assert np.allclose(vd.vis[:, 0], exp, rtol=1e-3, atol=1e-3)
    vn = VisibilitySimulator(add_noise=True, seed=3).run(sky, obs)
    assert vn.weight[0, 0] == pytest.approx(vd.weight[0, 0])


def test_visibility_simulator_agrees_with_iq_simulation():
    sky = SkyModel([PointSource("p", 0.0, 0.0, 2000.0)], center=FieldCenter("p", 120.0, 45.0))
    site = Site.from_toml(CFG / "sites" / "wilhelmsburg.toml")
    from interferosim.instrument import Antenna, Receiver
    dish = Antenna.from_toml(CFG / "antennas" / "dish_3m.toml")
    rec = Receiver.from_toml(CFG / "receivers" / "b210.toml")
    arr = Array(site, [Station.from_enu("A", site, 0, 0, 0, dish),
                       Station.from_enu("B", site, 200, 50, 0, dish)], receiver=rec)
    iq = IQSimulator(sky, arr, sample_rate_sps=2e4, seed=2)
    m = measure_visibilities(iq.simulate(20.0), iq, 20.0)
    vd_iq = VisData.from_measurements(m, iq.fc)
    assert abs(vd_iq.vis[0, 0] - 2000.0) < 5 * vd_iq.sigma[0, 0]


def test_transit_sun_two_horns():
    ra, dec = sun_radec_deg(mjd_from_iso("2026-11-01T11:00:00"))
    sun = SkyModel.from_json(CFG / "skies" / "sun_21cm.json")
    sun.center = FieldCenter("Sun", float(ra), float(dec))
    arr = Array.from_toml(CFG / "sites" / "two_horns_10m.toml")
    vd = TransitSimulator(arr, dec_deg=float(dec), ha_range_h=(-3, 3), t_dump_s=10.0, add_noise=False).run(sun)
    a = np.abs(vd.vis[:, 0])
    k = int(np.argmax(a))
    assert abs(vd.ha_deg[k]) < 1.0                       # maximum at the meridian
    assert a[k] == pytest.approx(abs(sun.visibility(vd.u[k], vd.v[k])), rel=0.01)
    assert a[0] < 0.2 * a[k]                             # drifted out of the horn beam (34 deg)
    ph = np.unwrap(np.angle(vd.vis[:, 0]))
    assert abs(ph[-1] - ph[0]) > 2 * np.pi               # fringes run through


# ---------------------------------------------------------------- visdata
def test_visdata_roundtrip_and_selection(tmp_path):
    sky = SkyModel.from_json(CFG / "skies" / "cyg_a.json")
    arr = Array.from_toml(CFG / "sites" / "aas6_sim_stations.toml")
    obs = Observation.from_sky(sky, arr, ha_range_h=(-1, 1), t_int_s=900)
    vd = VisibilitySimulator(seed=4).run(sky, obs)
    p = vd.save_npz(tmp_path / "v.npz")
    back = VisData.load_npz(p)
    assert np.allclose(back.vis, vd.vis) and back.header["OBJECT"] == vd.header["OBJECT"]
    assert vd.select_baseline("A-B").n_rows == vd.n_rows
    assert vd.average_channels().n_chan == 1
    with pytest.raises(ValueError):
        VisData(time_mjd=[0, 1], ha_deg=[0], u=[0, 1], v=[0, 1], w=[0, 1], baseline=[0, 0],
                vis=[1, 2], weight=[1, 1], freq_hz=[1e9])
