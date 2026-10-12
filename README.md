<!-- README.md -- Interferosim project overview -->
# Interferosim

An open-source simulator for radio interferometry — from a two-element
system up to small arrays — with interactive Jupyter notebooks for
teaching. A project of the University of Vienna.

Repository: <https://github.com/eschubi/interferosim> (private until the
licence questions are settled).

**Status:** 0.2.1 — sky models implemented (`interferosim.sky`, maps in
`interferosim.plotting`, notebook `notebooks/01_sky_models.ipynb`). The
further stages are being migrated from `aas6sim` and the teaching
notebooks step by step.

## Quick start

Create the project environment `.venv` (Python >= 3.11 required):

```
setup_venv.bat                 # Windows (double-click or command line)
./setup_venv.sh                # Linux / macOS
```

Then:

```
.venv\Scripts\activate          # Windows   (Linux/macOS: source .venv/bin/activate)
interferosim-sky configs/skies/aas6_demo.json --overview -o demo.png
jupyter lab notebooks          # kernel "Python (interferosim)"
```

## Processing chain

| Stage | Package module |
|---|---|
| Sky models (point, extended, HI 21 cm, OH masers, solar disc) | `interferosim.sky` |
| Antennas (selectable types), receiver (Ettus USRP B210), array | `interferosim.instrument` |
| Geometry, uv-coverage, visibilities, correlator | `interferosim.observation` |
| Gridding, weighting, dirty image and beam | `interferosim.imaging` |
| Högbom-CLEAN, Wiener filter, restoration | `interferosim.deconvolution` |
| FITS input/output | `interferosim.io` |

## Repository layout

```
src/interferosim/   Python package (core, independent of the user interface)
configs/            antenna types, receivers, sites (TOML)
notebooks/          one notebook per topic
tests/              pytest and nbval
binder/             environment for Binder / JupyterHub
docs/               LaTeX documentation (build: cd docs && latexmk)
```

## Installation (development, without the scripts)

```
python -m venv .venv
.venv\Scripts\python -m pip install -e ".[notebooks,test]"
.venv\Scripts\python -m pytest
```

Requires Python >= 3.11.

## API documentation

The source code documentation is generated with Doxygen
(<https://www.doxygen.nl>):

```
doxygen Doxyfile                    # -> docs/api/html/index.html
```

On Windows `build_api_docs.bat` does everything, including the
installation of Doxygen with winget.

## Licence

Creative Commons Attribution 4.0 International (CC BY 4.0), see `LICENSE`.

## Author

DI Erich Schubert
