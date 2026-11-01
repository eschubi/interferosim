# tests/test_import.py -- smoke test of the package skeleton
import importlib

import pytest

SUBPACKAGES = ["sky", "instrument", "observation", "imaging",
               "deconvolution", "io", "plotting", "widgets"]


def test_version():
    import interferosim
    assert interferosim.__version__ == "0.9.6"


@pytest.mark.parametrize("name", SUBPACKAGES)
def test_subpackage_imports(name):
    importlib.import_module(f"interferosim.{name}")


def test_constants():
    from interferosim import constants as c
    assert abs(c.wavelength(c.HI_21CM) - 0.2111) < 1e-4
    assert set(c.LINES) == {"HI", "OH1612", "OH1665", "OH1667", "OH1720"}


def test_configs_parse():
    import pathlib
    import tomllib
    root = pathlib.Path(__file__).resolve().parents[1] / "configs"
    files = sorted(root.rglob("*.toml"))
    assert files, "no configuration files found"
    for f in files:
        with open(f, "rb") as fh:
            tomllib.load(fh)
