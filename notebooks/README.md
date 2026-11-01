<!-- notebooks/README.md -->
# Notebooks

Teaching notebooks, one per topic. Planned notebooks and proposals are listed with their status; the migration of Appendix A (Wiener filter vs. CLEAN) and Appendix B (the Sun as a first target) becomes notebooks 07 and 08.

## Simulation chain

| Notebook | Content | Status |
|---|---|---|
| `00_general.ipynb` | project overview: goal of a VLB interferometer with an 11.5 km baseline (Wilhelmsburg–Michelbach), principle, key figures (resolution, delay, fringe rate, data rate), clock and timing requirements, sensitivity, road map | available |
| `01_sky_models.ipynb` | sky models, file format, maps, own skies, visibilities, catalogue of known regions with selector | available |
| `02_instrument.ipynb` | antenna types and primary beams, dish transfer function, B210 receiver and sensitivity, arrays and baselines, choice of the baseline length (Cyg A), selector | available |
| `03_signal_simulation.ipynb` | generating SDR receiver signals (IQ samples of the B210 at each station) based on the sky model under variable baseline parameters: signal model, raw fringes, measured vs. model visibilities over the baseline length, clock errors (TCXO/GPSDO), quantisation, SigMF output, selector | available |
| `04_observation.ipynb` | source tracking, geometric delay and fringe rate, uv-tracks by Earth rotation, visibilities over a track, transit (drift scan) of fixed antennas (Sun with two horns), asymmetry of a track, selector | available |
| `05_correlator.ipynb` | FX correlator on the simulated SDR signals: fractional delay, fringe stopping, delay/rate search and refinement, cross spectra, visibilities; station clock and LO offsets; selector | available |
| `06_imaging.ipynb` | gridding, natural/uniform/Wiener weighting, dirty beam and dirty image: OH masers with the project baseline and three stations, Cygnus A, selector | available |
| `07_deconvolution.ipynb` | Högbom CLEAN, clean beam and restoration, Wiener filter (migrated from Appendix A): maser group with the project baseline, convergence and loop gain, Wiener against CLEAN, three stations, limits for Cygnus A, selector | available |
| `08_sun_first_target.ipynb` | the Sun with two horns on one B210 (migrated from Appendix B): signal strength, baseline choice, fringe rate over the year, raw-signal checks (covered horns, crosstalk, channel phase), transit with the solar day, radio diameter, procedure, selector | available |
| `09_spectral_lines.ipynb` | spectral lines: channels and velocities, Doppler correction to the LSR, sky models with spectra, OH masers on the 11.5 km baseline (spectra towards a position, channel maps), effect of a missing Doppler correction, HI absorption against Cygnus A, planning of a line observation, selector | available |

## Theory

| Notebook | Content | Status |
|---|---|---|
| `A_fringes.ipynb` | shows and explains how fringes evolve: fringe pattern on the sky, fringe rate from the Earth's rotation, effect of the source size (fringe contrast), fringes in a long track (Cygnus A with zooms of 5–20 min, natural and structure fringe rate), selector | available (first version) |
| `B_fourier.ipynb` | explains the Fourier transformation and the correlation: transform pairs and convolution theorem, van Cittert–Zernike (numerical check against the analytic visibilities), one baseline = one Fourier component, the dish as low-pass filter, cross-correlation and cross spectrum (XF vs. FX, Wiener–Khinchin), the delay beam, selector | available |
| `C_sampling_iq.ipynb` | sampling and complex baseband: Nyquist and aliasing, IQ mixing, IQ imbalance and DC offset, band edges, quantisation (efficiency, input level, Van Vleck), what the B210 delivers, selector | available |
| `D_noise_sensitivity.ipynb` | radiometer equation, SEFD and antenna temperature, noise of a visibility (simulation check), S/N against integration time, bandwidth and baseline length for the project targets, coherence and clocks (random-walk phase, coherent vs. incoherent averaging), image noise, selector | available |
| `E_coordinates.ipynb` | coordinate frames: equatorial (RA/Dec, hour angle, sidereal time), horizontal (azimuth/elevation, circumpolar sources, parallactic angle), WGS84/ECEF/ENU and the baseline vector, uvw and the uv-ellipse (checked against the classical formulas), delay and w, precession J2000 → date, selector | available |

*available*: in the repository; *planned*: decided, not yet written;
*proposal*: suggestion for discussion. The former Appendices A (Wiener
filter vs. CLEAN) and B (the Sun) of the teaching material become the
numbered notebooks 07 and 08, so that the letters A, B, C, … are free for
the theory notebooks.
