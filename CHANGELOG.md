<!-- CHANGELOG.md -- Interferosim change log -->
# Changelog

All versions use three-part numbers (major.minor.patch); every release is
strictly higher than the previous one.

## [Unreleased]

### Added
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
  simulation); longitude and altitude of `wilhelmsburg.toml`.
- `binder/postBuild`; `repository-code` in `CITATION.cff`.

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
