# tests/test_solar.py -- tests of interferosim.observation.solar (0.8.1)
"""Sun flux, solar transit, solar-day transit simulation, disc fit."""

from pathlib import Path

import numpy as np
import pytest

from interferosim.instrument import Array
from interferosim.observation import (SolarTransitSimulator, TransitSimulator, disc_amplitude,
                                      fit_disc_diameter, gmst_deg, mjd_from_iso, solar_transit_mjd,
                                      sun_flux_jy, sun_radec_deg)
from interferosim.sky import FieldCenter, SkyModel

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs"
LAM = 299_792_458.0 / 1.420405751768e9


def _horns():
    return Array.from_toml(CFG / "sites" / "two_horns_10m.toml")


def _sun():
    return SkyModel.from_json(CFG / "skies" / "sun_21cm.json")


def test_sun_flux_matches_config():
    # T_b = 1e5 K, 34 arcmin -> the 476 200 Jy of sun_21cm.json (notebook B)
    assert sun_flux_jy(1e5, 34 / 60) == pytest.approx(476_200, rel=1e-4)


def test_solar_transit_is_local_noon():
    lon = _horns().reference.lon_deg
    t = solar_transit_mjd("2026-11-01", lon)
    ra, _ = sun_radec_deg(t)
    lha = (float(gmst_deg((t - 40587.0) * 86400.0)) + lon - float(ra) + 180) % 360 - 180
    assert lha == pytest.approx(0.0, abs=0.01)
    # equation of time on 1 November: apparent noon about 16 min before mean noon
    mean_noon = mjd_from_iso("2026-11-01") + 0.5 - lon / 360.0
    assert (mean_noon - t) * 1440 == pytest.approx(16.4, abs=0.5)


def test_solar_transit_simulator_amplitude_and_solar_rate():
    horns, sun = _horns(), _sun()
    a = SolarTransitSimulator(horns, "2026-11-01", (-2, 2), 5.0, add_noise=False).run(sun)
    # peak of |V| at transit: 476 kJy x |V(10 m)| of the limb-brightened disc
    q = 10.0 / LAM
    v10 = abs(sun.components[0].visibility(np.array([q]), np.array([0.0]))[0])
    assert np.abs(a.vis[:, 0]).max() == pytest.approx(v10, rel=0.01)
    k = int(np.argmax(np.abs(a.vis[:, 0])))
    assert abs(a.ha_deg[k]) < 0.1
    # same turns per hour angle as the sidereal model, but 0.27 % more time per turn
    t = solar_transit_mjd("2026-11-01", horns.reference.lon_deg)
    ra, dec = (float(x) for x in sun_radec_deg(t))
    sun.center = FieldCenter("Sun", ra, dec)
    b = TransitSimulator(horns, dec_deg=dec, ha_range_h=(-2, 2), t_dump_s=5.0, add_noise=False,
                         date="2026-11-01").run(sun)

    def rate(vd):
        ph = np.unwrap(np.angle(vd.vis[:, 0]))
        tt = (vd.time_mjd - vd.time_mjd[0]) * 86400.0
        return (ph[-1] - ph[0]) / (tt[-1] - tt[0])

    assert rate(b) / rate(a) == pytest.approx(1.00274, abs=3e-4)


def test_disc_fit_recovers_uniform_disc():
    q = np.linspace(5, 180, 25)
    amp = disc_amplitude(q, 4e5, 0.55)
    s0, d, es, ed = fit_disc_diameter(q, amp, p0=(3e5, 0.5))
    assert s0 == pytest.approx(4e5, rel=1e-6)
    assert d == pytest.approx(0.55, rel=1e-6)
    # first null at x = 3.8317
    assert disc_amplitude(3.8317 / (np.pi * np.radians(0.55)), 1.0, 0.55) < 1e-4
