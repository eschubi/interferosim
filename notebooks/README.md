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
| `04_observation.ipynb` | source tracking, geometric delay and fringe rate, uv-tracks by Earth rotation, visibilities over a track, transit (drift scan) of fixed antennas (Sun with two horns), selector | available |
| `05_correlator.ipynb` | FX correlator on the simulated SDR signals: fractional delay, fringe stopping, delay/rate search, cross spectra, visibilities | proposal |
| `06_imaging.ipynb` | gridding, natural/uniform/robust weighting, dirty beam and dirty image | proposal |
| `07_deconvolution.ipynb` | Högbom CLEAN vs. Wiener filter, restoration (migrated from Appendix A) | proposal |
| `08_sun_first_target.ipynb` | the Sun with two horns on one B210: transit, fringes, disc visibility (migrated from Appendix B) | proposal |
| `09_spectral_lines.ipynb` | multichannel spectra: HI 21 cm and the OH maser lines, Doppler/LSR velocity, line visibilities per channel | proposal |

## Theory

| Notebook | Content | Status |
|---|---|---|
| `A_fringes.ipynb` | shows and explains how fringes evolve (geometric delay, baseline projection, fringe rate, effect of source size) | planned |
| `B_fourier.ipynb` | explains the details of the Fourier transformation and the correlation (van Cittert–Zernike theorem, cross-correlation vs. cross spectrum, FX vs. XF) | planned |
| `C_sampling_iq.ipynb` | sampling and complex baseband: Nyquist, IQ mixing, aliasing, quantisation, what the B210 delivers | proposal |
| `D_noise_sensitivity.ipynb` | radiometer equation, SEFD, noise of a visibility, S/N against integration time, bandwidth and baseline length | proposal |
| `E_coordinates.ipynb` | coordinate frames: equatorial, horizontal, WGS84/ECEF, local east-north-up, uvw | proposal |

*available*: in the repository; *planned*: decided, not yet written;
*proposal*: suggestion for discussion. The former Appendices A (Wiener
filter vs. CLEAN) and B (the Sun) of the teaching material become the
numbered notebooks 07 and 08, so that the letters A, B, C, … are free for
the theory notebooks.
