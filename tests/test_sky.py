# tests/test_sky.py -- tests of interferosim.sky and interferosim.plotting.sky
import json
import math
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from interferosim.sky import (AM, AS, DiscSource, FieldCenter, GaussianSource,  # noqa: E402
                              PointSource, SkyModel, component_from_dict)

ROOT = Path(__file__).resolve().parents[1]

# Excerpt of AAS6 json/interfero/sim_sky.json (German note key "_hinweis")
AAS6_SIM_SKY = {
    "_hinweis": "Simulationshimmel des AAS6-Interferometers",
    "center": {"label": "W3(OH)", "ra_deg": 36.7661, "dec_deg": 61.8742},
    "components": [
        {"name": "HII-Region A", "dx_deg": -2.0, "dy_deg": 1.2, "flux_jy": 900.0,
         "kind": "gauss", "fwhm_x_deg": 0.7, "fwhm_y_deg": 0.7, "pa_deg": 0.0},
        {"name": "Maser M1", "dx_deg": 0.3, "dy_deg": 0.45, "flux_jy": 70.0,
         "kind": "point", "fwhm_x_deg": 0.0, "fwhm_y_deg": 0.0, "pa_deg": 0.0},
    ],
}


# ---------------------------------------------------------------- format
def test_read_aas6_format(tmp_path):
    p = tmp_path / "sim_sky.json"
    p.write_text(json.dumps(AAS6_SIM_SKY), encoding="utf-8")
    sky = SkyModel.from_json(p)
    assert sky.center.label == "W3(OH)"
    assert sky.note.startswith("Simulationshimmel")
    assert isinstance(sky["HII-Region A"], GaussianSource)
    assert isinstance(sky["Maser M1"], PointSource)
    assert sky.total_flux_jy == pytest.approx(970.0)


def test_json_roundtrip_keeps_aas6_keys(tmp_path):
    sky = SkyModel.demo()
    sky.add(DiscSource("disc", 1.0, -1.0, 100.0, diameter_deg=0.5, limb=0.5))
    p = sky.to_json(tmp_path / "s.json")
    raw = json.loads(p.read_text(encoding="utf-8"))
    for c in raw["components"]:
        assert {"name", "dx_deg", "dy_deg", "flux_jy", "kind"} <= set(c)
    gauss = [c for c in raw["components"] if c["kind"] == "gauss"][0]
    assert {"fwhm_x_deg", "fwhm_y_deg", "pa_deg"} <= set(gauss)
    back = SkyModel.from_json(p)
    assert [c.to_dict() for c in back] == [c.to_dict() for c in sky]
    assert back.center == sky.center


def test_unknown_kind_and_missing_field():
    with pytest.raises(ValueError):
        component_from_dict({"name": "x", "dx_deg": 0, "dy_deg": 0, "flux_jy": 1,
                             "kind": "ring"})
    with pytest.raises(ValueError):
        component_from_dict({"name": "x", "dx_deg": 0, "flux_jy": 1})


def test_demo_matches_aas6sim():
    sky = SkyModel.demo()
    assert len(sky) == 9
    assert sky.total_flux_jy == pytest.approx(7695.0)
    m2 = sky["Maser M2"]
    assert m2.dx_deg == pytest.approx(0.30 + 5.2 * AS)
    assert m2.dy_deg == pytest.approx(0.45 - 1.5 * AS)
    assert sky["Compact source K1"].fwhm_x_deg == pytest.approx(2.1 * AM)


@pytest.mark.parametrize("name", ["aas6_demo", "w3oh_masers", "sun_21cm", "cas_a",
                                  "cyg_a", "tau_a", "vir_a", "orion_a"])
def test_example_skies_load(name):
    sky = SkyModel.from_json(ROOT / "configs" / "skies" / f"{name}.json")
    assert len(sky) > 0 and sky.total_flux_jy > 0


# ---------------------------------------------------------------- render
def test_render_conserves_flux():
    sky = SkyModel([GaussianSource("g", 0.0, 0.0, 100.0, 1.0, 0.5, pa_deg=30.0),
                    PointSource("p", 1.0, 1.0, 10.0),
                    DiscSource("d", -1.5, 0.0, 50.0, diameter_deg=0.8, limb=0.5)])
    img, x, y, pix = sky.render_map(fov_deg=6.0, npix=600, center_deg=(0.0, 0.0))
    assert img.sum() == pytest.approx(160.0, rel=2e-3)


def test_point_source_pixel_position():
    sky = SkyModel([PointSource("p", 0.5, -0.25, 7.0)])
    img, x, y, _ = sky.render_map(fov_deg=2.0, npix=200, center_deg=(0.0, 0.0))
    j, i = np.unravel_index(np.argmax(img), img.shape)
    assert x[i] == pytest.approx(0.5, abs=0.01)
    assert y[j] == pytest.approx(-0.25, abs=0.01)
    assert img.max() == pytest.approx(7.0)


def test_map_contains_offset_components():
    sky = SkyModel.from_json(ROOT / "configs" / "skies" / "w3oh_masers.json")
    img, x, y, _ = sky.render_map()
    assert img.sum() == pytest.approx(sky.total_flux_jy)


# ---------------------------------------------------------------- visibility
def test_visibility_at_origin_is_total_flux():
    sky = SkyModel.demo()
    sky.add(DiscSource("disc", 0.0, 0.0, 100.0, diameter_deg=0.5, limb=0.5))
    assert sky.visibility(0.0, 0.0).real == pytest.approx(sky.total_flux_jy, rel=1e-6)


def test_gaussian_visibility_formula():
    theta = 0.1                                   # FWHM in degrees
    g = GaussianSource("g", 0.0, 0.0, 1.0, theta, theta)
    q = 50.0
    expected = math.exp(-math.pi**2 * math.radians(theta)**2 * q**2 / (4 * math.log(2)))
    assert abs(g.visibility(q, 0.0)) == pytest.approx(expected, rel=1e-9)


def test_gaussian_visibility_matches_fft_of_render():
    g = GaussianSource("g", 0.0, 0.0, 1.0, 0.3, 0.15, pa_deg=40.0)
    sky = SkyModel([g])
    n, pix = 512, 0.01
    img, x, y, _ = sky.render_map(fov_deg=n * pix, pix_deg=pix, center_deg=(0.0, 0.0))
    F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))
    du = 1.0 / (n * math.radians(pix))
    k = 20
    c = n // 2
    for (ju, jv) in [(k, 0), (0, k), (k, k)]:
        u, v = ju * du, jv * du
        # FFT convention exp(-2 pi i (u l + v m)); array index [v, u]
        assert abs(F[c + jv, c + ju]) == pytest.approx(abs(g.visibility(u, v)), abs=2e-3)


def test_uniform_disc_first_null():
    d = DiscSource("d", 0.0, 0.0, 1.0, diameter_deg=34.0 / 60.0)
    theta = math.radians(34.0 / 60.0)
    q_null = 1.2197 / theta
    assert abs(d.visibility(q_null, 0.0)) < 1e-3
    lb = DiscSource("d", 0.0, 0.0, 1.0, diameter_deg=34.0 / 60.0, limb=0.0)
    lb_num = DiscSource("d", 0.0, 0.0, 1.0, diameter_deg=34.0 / 60.0, limb=1e-12)
    for q in (20.0, 60.0, 100.0):     # Hankel path agrees with 2 J1(x)/x
        assert lb_num.visibility(q, 0.0).real == pytest.approx(
            lb.visibility(q, 0.0).real, abs=2e-4)


def test_point_source_phase():
    p = PointSource("p", 0.01, 0.0, 2.0)
    u = 1.0 / math.radians(0.01) / 4.0          # quarter turn
    assert np.angle(p.visibility(u, 0.0)) == pytest.approx(-math.pi / 2)


# ---------------------------------------------------------------- plotting
def test_plot_functions_run():
    from interferosim.plotting import plot_sky_map, plot_sky_overview
    sky = SkyModel.demo()
    for scale in ("log", "linear", "asinh"):
        fig, ax, img = plot_sky_map(sky, npix=128, scale=scale)
        assert img.shape == (128, 128)
    fig, axes = plot_sky_overview(sky, npix=128)
    assert len(axes) == 2


def test_plot_gallery_runs():
    from interferosim.plotting import plot_sky_gallery
    skies = [SkyModel.from_json(f) for f in sorted((ROOT / "configs" / "skies").glob("*.json"))]
    fig, axes = plot_sky_gallery(skies, npix=64, log_names={"aas6_demo"})
    assert axes.size >= len(skies)


def test_cli_writes_png(tmp_path):
    from interferosim.sky.cli import main
    out = tmp_path / "demo.png"
    assert main(["--demo", "--overview", "--npix", "128", "-o", str(out)]) == 0
    assert out.stat().st_size > 1000


# Perley & Butler (2017) at 1.42 GHz; Cas A and Tau A corrected to epoch 2026
@pytest.mark.parametrize("name, flux", [("cas_a", 1636.0), ("cyg_a", 1556.0),
                                        ("tau_a", 806.0), ("vir_a", 210.0)])
def test_a_team_total_flux(name, flux):
    sky = SkyModel.from_json(ROOT / "configs" / "skies" / f"{name}.json")
    assert sky.total_flux_jy == pytest.approx(flux, rel=1e-3)
    assert sky.center.ra_deg is not None and sky.note
