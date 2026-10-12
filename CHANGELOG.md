<!-- CHANGELOG.md -- Interferosim change log -->
# Changelog

All versions use three-part numbers (major.minor.patch); every release is
strictly higher than the previous one.

## [Unreleased]

## [0.3.1] — 2026-10-12

### Added
- Documentation, chapter "Stations and Array": section "Choice of the
  baseline length" with the example Cyg A (|V| and noise over the baseline,
  three regimes, options for the long baseline), figure `cyga_baseline.png`.
- Notebook 02, section 3: the same analysis as code (figure and table).

## [0.3.0] — 2026-10-12

### Added
- Package `interferosim.instrument`:
  - `Antenna` (abstract) with `Horn`, `Yagi`, `ParabolicDish`: gain,
    effective area, HPBW, primary beam, `jy_to_k`, `from_toml`; dish with
    `u_max`, `u_cut_eff`, `otf` (from aas6sim `model.Instrument`).
  - `Receiver` (B210 profile): `sefd_jy`, `single_dish_sigma_jy`,
    `baseline_sigma_jy`.
  - `Site`, `Station`, `Array`: WGS84 <-> ECEF (from AAS6
    `interfero/geometry.py`), local east/north/up frame, baselines;
    `Array.from_toml` reads site, antennas and receiver named in the file.
- `interferosim.plotting`: `plot_beams`, `plot_array`.
- Notebook `02_instrument.ipynb` (antennas, dish transfer function,
  sensitivity of both arrays for the sky models, arrays, selector).
- `configs/sites/two_horns_10m.toml` (notebook B); `aas6_sim_stations.toml`
  names site, antenna and receiver.
- Documentation: chapters "Antennas" (implementation), "Receiver",
  "Stations and Array" (new), "Interactive Notebooks"; data formats and
  UML diagram updated.

### Changed
- `configs/antennas/yagi.toml`: typical example values (13 elements,
  0.85 m boom, 15 dBi) instead of `nan` (decision 2026-10-12).

## [0.2.2] — 2026-10-12

### Added
- Known sky regions at 1420 MHz in `configs/skies/`: `cas_a`, `cyg_a`,
  `tau_a`, `vir_a` (flux densities from Perley & Butler 2017, Cas A/Tau A
  corrected for secular fading to epoch 2026) and `orion_a` (approximate).
- Notebook 01, section 6: selector (sky file, scale, view incl. visibility
  amplitude over the baseline) and a static gallery of all sky files.
- `interferosim.plotting.plot_sky_gallery()`: several sky models side by
  side; used in notebook 01 and for the documentation figure
  `docs/images/sky_gallery.png`.
- README: section "Example skies".
- `Doxyfile`: API documentation of the package with Doxygen
  (HTML and LaTeX into `docs/api/`, excluded from version control).
- Documentation: chapters "Version Control", "API Documentation",
  "Tests and Validation" (running the tests), git workflow diagram.
- Documentation: requirements with identifiers, chapter "Data Formats",
  antenna types (Yagi; proposals helix, LPDA, dipole, patch), basic design
  0.2 (package `signal` for the closed simulation with SDR signal,
  spectra per component, `TransitSimulator`, module `visdata`).
- Configuration: `configs/antennas/yagi.toml` (template),
  `configs/sites/aas6_sim_stations.toml` (stations A/B of the AAS6
  simulation); `wilhelmsburg.toml` with the confirmed position of station A
  (48.088056°, 15.756111°, 470 m).
- `binder/postBuild`; `repository-code` in `CITATION.cff`.
- `build_api_docs.bat` / `build_api_docs.ps1`: install Doxygen (winget) if
  needed, generate and open the API documentation.

## [0.2.1] — 2026-10-12

### Changed
- Demo sky (`SkyModel.demo()`, `configs/skies/aas6_demo.json`): component
  names in English (positions and fluxes unchanged; the AAS6 names are
  listed in the docstring of `SkyModel.demo`).

### Added
- Project virtual environment scripts: `setup_venv.ps1` / `setup_venv.bat`
  (Windows) and `setup_venv.sh` (Linux/macOS); they create `.venv`, install
  the package with the extras `notebooks,test`, register the Jupyter kernel
  "Python (interferosim)" and run the tests.
- `ipykernel` added to the extra `notebooks`.

## [0.2.0] — 2026-10-12

### Added
- Package `interferosim.sky`: `SkyComponent` (abstract), `PointSource`,
  `GaussianSource` (rendering and visibility from aas6sim 1.3.1),
  `DiscSource` (uniform / limb-brightened disc, notebook B),
  `FieldCenter`, `SkyModel` with `render_map()`, `visibility()`,
  `demo()` and JSON input/output in the format of AAS6
  `json/interfero/sim_sky.json` (new component kind `disc`).
- `interferosim.plotting`: `plot_sky_map()`, `plot_sky_overview()`.
- Command line tool `interferosim-sky`.
- Example skies `configs/skies/` (aas6_demo, w3oh_masers, sun_21cm).
- Notebook `notebooks/01_sky_models.ipynb`.
- Tests `tests/test_sky.py` (format, flux conservation, visibilities
  against FFT and analytic formulas, plots, CLI).
- Documentation: chapter "Sky Models", UML class diagram updated.

## [0.1.0] — 2026-10-12

### Added
- Repository skeleton: `pyproject.toml`, `LICENSE` (CC BY 4.0),
  `CITATION.cff`, `README.md`, `.gitignore`, `binder/environment.yml`.
- Package skeleton `src/interferosim/` with sub-packages `sky`,
  `instrument`, `observation`, `imaging`, `deconvolution`, `io`,
  `plotting`, `widgets`.
- Configuration files: antenna types (horn, 3 m dish), receiver (B210),
  sites.
- Smoke test `tests/test_import.py`.
- Documentation moved to an English, self-contained LaTeX project
  (`docs/interferosim.tex`, local style copy, project-specific list of
  abbreviations, `latexmkrc`).
