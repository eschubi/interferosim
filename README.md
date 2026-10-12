<!-- README.md -- Interferosim project overview -->
# Interferosim

An open-source simulator for radio interferometry — from a two-element
system up to small arrays — with interactive Jupyter notebooks for
teaching. A project of the University of Vienna.

Repository: <https://github.com/eschubi/interferosim> (private until the
licence questions are settled).

**Status:** 0.3.1 — sky models (`interferosim.sky`, notebook
`01_sky_models.ipynb` with a selector for all example skies) and
instrument (`interferosim.instrument`: horn, parabolic dish and Yagi,
Ettus USRP B210 receiver, sites, stations and arrays; notebook
`02_instrument.ipynb`) implemented. The
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

## Example skies

`configs/skies/` contains sky models in the JSON format of AAS6
(`json/interfero/sim_sky.json`, plus the component kind `disc`):

| File | Content | S (1.42 GHz) |
|---|---|---|
| `aas6_demo.json` | demo field of `aas6sim` around W3(OH) | 7695 Jy |
| `w3oh_masers.json` | the three OH masers of the demo field | 135 Jy |
| `sun_21cm.json` | quiet Sun, limb-brightened disc 34′ | 476 kJy |
| `cas_a.json` | Cassiopeia A, supernova remnant | 1636 Jy |
| `cyg_a.json` | Cygnus A, radio galaxy (lobes, hotspots, core) | 1556 Jy |
| `tau_a.json` | Taurus A (Crab Nebula) | 806 Jy |
| `vir_a.json` | Virgo A (M87): core, jet, inner lobes, halo | 210 Jy |
| `orion_a.json` | Orion A (M42), approximate model | ~400 Jy |

Flux densities of Cas A, Cyg A, Tau A and Vir A follow Perley & Butler
(2017, ApJS 230, 7); Cas A and Tau A are corrected for their secular fading
to epoch 2026. Each file documents sources and simplifications in `_note`.
In notebook 01, section 6, a selector shows each sky as map, overview or
visibility amplitude over the baseline;
`interferosim.plotting.plot_sky_gallery()` draws all of them side by side.

![Example skies](docs/images/sky_gallery.png)

## Instrument

```python
from interferosim.instrument import Antenna, Receiver, Array
dish = Antenna.from_toml("configs/antennas/dish_3m.toml")
b210 = Receiver.from_toml("configs/receivers/b210.toml")
arr  = Array.from_toml("configs/sites/aas6_sim_stations.toml")   # site, dishes, B210
b210.sefd_jy(dish), arr.baselines_enu(), b210.baseline_sigma_jy(dish, dish, 10.0)
```

Antennas: `horn_15dbi`, `dish_3m`, `yagi` (typical example values);
arrays: `aas6_sim_stations` (two 3 m dishes, 11.1 km east-west) and
`two_horns_10m` (notebook B).

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
