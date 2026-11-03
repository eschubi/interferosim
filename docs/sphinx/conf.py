# docs/sphinx/conf.py -- Sphinx configuration of the Interferosim documentation
# Build:  build_sphinx_docs.bat   or   python -m sphinx -b html docs/sphinx docs/sphinx/_build/html
# Version 1.0.0 (2026-10-18): first version.

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import interferosim  # noqa: E402

# -- Project ---------------------------------------------------------------
project = "Interferosim"
author = "Erich Schubert"
copyright = "2026, Erich Schubert (MIT License)"
release = interferosim.__version__
version = ".".join(release.split(".")[:2])

# -- General ---------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "sphinx_copybutton",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build"]
language = "en"

# -- autodoc / napoleon ----------------------------------------------------
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autoclass_content = "both"

napoleon_numpy_docstring = True
napoleon_google_docstring = False
# Section headings used in the module docstrings of interferosim
napoleon_custom_sections = [
    ("Origin", "notes_style"),
    ("Conventions", "notes_style"),
    ("Formulas", "notes_style"),
    ("File format", "notes_style"),
    ('Relations (chapter "Radio Astronomical Fundamentals")', "notes_style"),
    ("Configuration (``configs/sites/*.toml``)", "notes_style"),
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

myst_enable_extensions = ["dollarmath"]

# -- HTML ------------------------------------------------------------------
html_theme = "furo"
html_title = f"Interferosim {release}"
html_static_path = []
