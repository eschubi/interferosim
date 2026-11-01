<!-- CHANGELOG.md -- Interferosim change log -->
# Changelog

All versions use three-part numbers (major.minor.patch); every release is
strictly higher than the previous one.

## [Unreleased]

## [0.9.5] — 2026-11-01

### Changed
- Notebook `D_noise_sensitivity.ipynb`, section 7 "Summary": a formula table
  (step, quantity, formula, dependence, package function, section) and a
  table of project values computed in the notebook (conditions stated per
  row, 12- and 2-bit times, coherence, incoherent detection, image noise)
  instead of the former static table with rounded values.

## [0.9.4] — 2026-11-01

### Added
- `interferosim.observation.geometry`: `lst_deg`, `azimuth_elevation_deg`,
  `parallactic_angle_deg`, `precess_deg` (IAU 1976 precession, within 0.1″
  of astropy).
- Theory notebook `E_coordinates.ipynb`; with it all planned notebooks
  (00–09, A–E) are available.
- Documentation: subsection "Horizontal coordinates and precession" in
  chapter "Observation" (one figure); test.

## [0.9.3] — 2026-11-01

### Added
- Module `interferosim.instrument.sensitivity`: `antenna_temperature_k`,
  `radiometer_sigma_k`, `visibility_sigma_jy`, `integration_time_s`,
  `image_noise_jy`, `coherence_factor`, `phase_rms_white_fm_rad`,
  `coherence_time_s`, `random_walk_coherence`, `incoherent_snr`.
- Theory notebook `D_noise_sensitivity.ipynb`.
- Documentation: section "Sensitivity" in chapter "Receiver" (one figure);
  test.

## [0.9.2] — 2026-11-01

### Added
- Module `interferosim.signal.baseband`: `iq_imbalance`,
  `image_rejection_db`, `quantizer_levels`, `quantization_efficiency`,
  `van_vleck_1bit`.
- Theory notebook `C_sampling_iq.ipynb`.
- Doxygen: `docs/style/doxygen_unicode.sty` (Unicode characters of the
  docstrings for the PDF reference manual), `LATEX_EXTRA_STYLESHEET` in the
  `Doxyfile`; the reference manual builds without errors (321 pages for 0.9.1).
- Documentation: section "Baseband and quantisation" in chapter "Signal
  Simulation" (two figures); tests.

## [0.9.1] — 2026-11-01

### Added
- `interferosim.signal`: `lag_correlation` (cross-correlation function over
  the lags) and `xf_cross_spectrum` (XF correlator, Wiener–Khinchin).
- Theory notebook `B_fourier.ipynb` (transform pairs, convolution theorem,
  van Cittert–Zernike with numerical check, one baseline = one Fourier
  component, aperture autocorrelation, XF vs. FX, delay beam, selector).
- Documentation: section "XF and FX" in chapter "Correlator" (two figures);
  test.

## [0.9.0] — 2026-11-01

### Added
- Module `interferosim.sky.spectra`: `FlatSpectrum`, `PowerLawSpectrum`,
  `GaussianLineSpectrum`, `AbsorptionSpectrum`, `spectrum_from_dict`,
  `radio_velocity_kms`, `frequency_from_velocity`.
- `SkyModel.spectra` (JSON key `"spectrum"` per component), `set_spectrum`,
  `flux_at`, `at_frequency`, `spectrum_jy`, `visibility_spectrum`.
- Module `interferosim.observation.doppler`: `observer_velocity_kms`,
  `v_lsr_correction_kms` (checked against astropy within 0.02 km/s),
  `topocentric_frequency_hz`, `lsr_velocity_from_topocentric`,
  `regrid_to_lsr`.
- `interferosim.observation.SpectralVisibilitySimulator` (multichannel
  visibilities with spectra and Doppler shift).
- `Imager.channel_images` (image cube).
- Sky models `w3oh_masers_lines.json`, `cyg_a_hi.json`.
- Notebook `09_spectral_lines.ipynb`.
- Documentation: chapter "Spectral Lines" (five figures); tests.

## [0.8.1] — 2026-10-31

### Added
- Module `interferosim.observation.solar`: `sun_flux_jy`,
  `solar_transit_mjd` (local apparent noon), `SolarTransitSimulator`
  (drift scan with the position of the Sun recomputed per dump, solar-day
  rate), `disc_amplitude`, `fit_disc_diameter`.
- Notebook `08_sun_first_target.ipynb` (Appendix B: signal strength,
  baseline choice, fringe rate over the year, raw-signal checks with the IQ
  simulation, transit, radio diameter, procedure, selector).
- Documentation: section "The Sun as the first target" in chapter
  "Observation" (five figures); tests.

## [0.8.0] — 2026-10-31

### Added
- Package `interferosim.deconvolution`: interface `Deconvolver`
  (`deconvolve`, alias `run`), `HogbomClean` (loop gain, iteration limit,
  absolute and relative threshold, CLEAN window, positivity; convergence
  history), `CleanBeam` (Gaussian fitted to the main lobe, kernel,
  convolution), `DeconvResult` (model, residual, restored image,
  components, flux within a circle), `WienerFilter` (W = S*/(|S|^2+NSR),
  optional target resolution), `reference_image` (sky model on the image
  grid, optionally attenuated by the primary beam and convolved with the
  clean beam), `rms_error`.
- `interferosim.plotting`: `plot_clean_history`.
- Notebook `07_deconvolution.ipynb` (principle, CLEAN on the maser group,
  convergence and loop gain, Wiener against CLEAN, three stations, limits
  for Cygnus A, selector).
- Documentation: chapter "Deconvolution" (implementation, results, five
  figures); Sphinx API page; tests.

## [0.7.0] — 2026-10-31

### Added
- Package `interferosim.imaging`: `Gridder` (nearest-cell gridding with
  Hermitian conjugates, multi-frequency synthesis, phase-centre shift),
  `NaturalWeighting`, `UniformWeighting`, `WienerWeighting` (notebook A),
  `Imager` (dirty image and dirty beam in Jy/beam), `Image` (peak, value at
  a position, rms, main-lobe FWHM).
- `interferosim.plotting`: `plot_image`, `sky_marks`.
- Notebook `06_imaging.ipynb` (principle, OH masers with the project
  baseline, weighting comparison, three stations, Cygnus A, selector).
- Documentation: chapter "Imaging" (implementation, results, four figures);
  Sphinx API page; tests.

### Changed
- Road map: imaging in 0.7, deconvolution in 0.8, spectral lines in 0.9.

## [0.6.0] — 2026-10-31

### Added
- FX correlator `interferosim.signal.FXCorrelator`: model correction,
  fringe search per station (`search`, with refinement), cross spectra per
  integration with weights from the auto spectra, result as `VisData` with
  a frequency axis.
- Correlation kernels taken over from AAS6 `modules/interfero/correlate.py`
  1.0.0: `integer_align`, `fringe_stop`, `Accumulator`, `cross_spectrum`,
  `FringeSolution`, `delay_rate_search`, `apply_fringe_solution`,
  `stokes_from_products`; new `refine_fringe_solution` (phase-slope and
  phase-drift fits remove the sub-sample bias of the search).
- Notebook `05_correlator.ipynb` (Sun with two dishes 15 m apart, each
  with its own B210; station B with clock and LO offset: cross and lag
  spectrum, search surface, refinement, cross spectra after the solution,
  channels and noise, selector).
- Documentation: chapter "Correlator" with three figures; Sphinx API page;
  tests for the kernels and the correlator.

### Changed
- Road map: correlator in version 0.6, spectral lines in 0.7, imaging and
  deconvolution in 0.8.

## [0.5.1] — 2026-10-31

### Added
- Theory notebook `A_fringes.ipynb` (first version): fringe pattern on the
  sky, fringe rate from the Earth's rotation, fringe contrast and source size
  (Sun with two horns), fringes in a long track -- Cygnus A on 11.5 km with
  zooms of 5, 10, 15 and 20 min and the natural and structure fringe rates;
  selector.
- Notebook 04: section "Why the track is not symmetric about transit" with
  the counter-check (jet axis turned east-west); documentation chapter
  "Observation": subsection "Asymmetry of a track", figure
  `obs_track_symmetry.png`.

## [0.5.0] — 2026-10-31

### Added
- `interferosim.observation`: `Observation` (hour angles, times of a
  transit on a date, elevation per station, visibility mask, uvw, delay,
  fringe rate), `UVCoverage`, `VisibilitySimulator` (analytic visibilities
  with radiometer noise, tracking antennas), `TransitSimulator` (fixed
  antennas, drift scan without fringe stopping, notebook B); module `sun`
  (Sun position, low precision; `mjd_from_iso`).
- Module `interferosim.visdata`: `VisData` (rows per integration and
  baseline, columns as aas6sim plus baseline index and frequency axis;
  `select_baseline`, `channel`, `average_channels`, `concatenate`,
  `save_npz`/`load_npz`, `from_measurements` for the results of notebook 03).
- `interferosim.plotting`: `plot_uv_coverage`, `plot_elevation`, `plot_visdata`.
- Notebook `04_observation.ipynb` (sources over the day, delay and fringe
  rate, uv-coverage for several declinations and a hypothetical
  three-station array, visibilities over a track, Sun transit with two
  horns, selector).
- Documentation: chapter "Observation" (observation and uv-coverage,
  simulated visibilities, transit), VisData in "Data Formats", UML and
  design text updated, figures; Sphinx API pages; tests
  `tests/test_observation.py`.

## [0.4.1] — 2026-10-31

### Changed
- Notebook 03, section 2: the titles of the IQ figure say which part of the
  recording is shown (first 2 ms of 10 s, all samples for the histogram,
  number of blocks for the spectrum); new slider to move a window of
  selectable length through the whole recording of both stations, with the
  correlation coefficient of the window.

## [0.4.0] — 2026-10-31

### Added
- Package `interferosim.signal`: `IQSimulator` (IQ samples of all stations
  from the sky model: covariance from the analytic component visibilities,
  phase-centre delay and fringe phase per sample, primary beams, receiver
  noise), `IQBlock` (with SigMF output), `StationErrors` (LO offset and
  phase, clock offset; `from_receiver` within the data-sheet limits),
  `quantize`; `correlate`: `fractional_delay` (from AAS6), `correct_geometry`,
  `cross_correlate`, `measure_visibilities`.
- `interferosim.observation.geometry`: GMST, Greenwich hour angle, sky basis
  and source direction in ECEF, geometric delay and rate, uvw, elevation.
- Notebook `03_signal_simulation.ipynb` (signal model, IQ statistics and
  spectrum, raw fringes, measured vs. model visibilities over the baseline
  length, clock errors, quantisation, SigMF, selector).
- Documentation: chapter "Signal Simulation", section "Geometry" in
  "Observation", figure `signal_check.png`; Sphinx API pages for
  `observation` and `signal`; tests `tests/test_signal.py`.

### Changed
- Road map: the signal simulation (package `signal`) is brought forward to
  version 0.4; observation (uv-coverage, transit) and `visdata` follow.

## [0.3.2] — 2026-10-31

### Added
- B210 device profile: `Receiver` with data-sheet fields (RF range, ADC,
  maximum sample rate and bandwidth for one and two channels, noise figure,
  gain range, IIP3, SFDR, reference accuracy TCXO/GPSDO, PPS accuracy, host
  interface, FPGA) and operating settings `center_freq_hz`, `rx_gain_db`,
  `sample_format`, `time_source`, `external_lna`; derived quantities
  `data_rate_bytes_per_s`, `receiver_noise_temperature_k`,
  `relative_frequency_error`, `lo_offset_hz`, `phase_coherence_time_s`;
  `validate()`; profiles with `base = "<name>"`.
- `configs/receivers/b210.toml` 0.2.0 (data-sheet values), new profiles
  `b210_gpsdo.toml` (one VLBI station with GPSDO) and `b210_aas6_sim.toml`
  (settings of the AAS6 simulation).
- `docs/images/b210_blockdiagram.tex` (TikZ, standalone) with PDF/PNG:
  receive path of the B210; notebook 02 section "B210 device profile";
  documentation chapter "Receiver" extended, data formats, bibliography.
- Hello-world example `interferosim.hello` / `interferosim-hello`
  (E. Schubert).
- `docs/images/array_geometry.tex` (TikZ, standalone class) and its PDF:
  array of J antennas with geometric delay, primary beams, field of view,
  baseline and receiver chains g_k -> x_k(t); PNG rendering embedded in
  notebook 00 (section 1) with an explanation.
- Notebook `00_general.ipynb`: project overview with the goal of a VLB
  interferometer on the 11.5 km baseline Wilhelmsburg-Michelbach (principle,
  VLBI concept, key figures, clock/timing requirements with B210 data-sheet
  values, sensitivity, road map).
- Sphinx documentation in `docs/sphinx/` (Furo theme): the radio
  interferometer in brief, installation, quick start, sky models,
  instrument, notebooks, API reference from the docstrings, change log;
  `build_sphinx_docs.bat`; optional dependency group `docs`.
- README, `notebooks/README.md`, Sphinx page "Notebooks": notebook plan --
  simulation chain 01-09 (03 signal simulation planned) and theory
  notebooks A-E (A fringes, B Fourier transformation and correlation
  planned), further notebooks as proposals.

### Removed
- `notebooks/01_sky_models.local-backup.ipynb.bak` (local backup copy,
  committed by mistake); `*.bak` added to `.gitignore`.

### Changed
- Notebooks 00-02: first cell with `id: notebooks/<name>.ipynb` and the short
  description of the notebook (01: the machine-specific localhost link
  replaced by the repository path); author line below the title of each
  notebook (as added by E. Schubert in notebook 00).

### Fixed
- Notebooks 00-02 name the kernel `interferosim` ("Python (interferosim)",
  installed by `setup_venv`) instead of the generic `python3`, so Jupyter
  starts them in the project environment; `binder/postBuild` registers the
  kernel.
- Docstrings of `DiscSource` and `otf_circular`: valid reStructuredText.

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
