.. docs/sphinx/notebooks.rst -- teaching notebooks

Notebooks
=========

The teaching notebooks in ``notebooks/`` use the package and add
interactive selectors (ipywidgets). Start them with

.. code-block:: text

   jupyter lab notebooks          # kernel "Python (interferosim)"

Simulation chain
----------------

.. list-table::
   :header-rows: 1
   :widths: 28 60 12

   * - Notebook
     - Content
     - Status
   * - ``01_sky_models.ipynb``
     - sky models, file format, maps, own skies, visibilities, catalogue of known regions with selector
     - available
   * - ``02_instrument.ipynb``
     - antenna types and primary beams, dish transfer function, B210 receiver and sensitivity, arrays and baselines, choice of the baseline length (Cyg A), selector
     - available
   * - ``03_signal_simulation.ipynb``
     - generating SDR receiver signals (IQ samples of both B210 channels) based on the sky model under variable baseline parameters
     - planned
   * - ``04_observation.ipynb``
     - source tracking, geometric delay and fringe rate, uv-tracks by Earth rotation, transit (drift scan) of fixed antennas
     - proposal
   * - ``05_correlator.ipynb``
     - FX correlator on the simulated SDR signals: fractional delay, fringe stopping, delay/rate search, cross spectra, visibilities
     - proposal
   * - ``06_imaging.ipynb``
     - gridding, natural/uniform/robust weighting, dirty beam and dirty image
     - proposal
   * - ``07_deconvolution.ipynb``
     - Högbom CLEAN vs. Wiener filter, restoration (migrated from Appendix A)
     - proposal
   * - ``08_sun_first_target.ipynb``
     - the Sun with two horns on one B210: transit, fringes, disc visibility (migrated from Appendix B)
     - proposal
   * - ``09_spectral_lines.ipynb``
     - multichannel spectra: HI 21 cm and the OH maser lines, Doppler/LSR velocity, line visibilities per channel
     - proposal

Theory
------

.. list-table::
   :header-rows: 1
   :widths: 28 60 12

   * - Notebook
     - Content
     - Status
   * - ``A_fringes.ipynb``
     - shows and explains how fringes evolve (geometric delay, baseline projection, fringe rate, effect of source size)
     - planned
   * - ``B_fourier.ipynb``
     - explains the details of the Fourier transformation and the correlation (van Cittert–Zernike theorem, cross-correlation vs. cross spectrum, FX vs. XF)
     - planned
   * - ``C_sampling_iq.ipynb``
     - sampling and complex baseband: Nyquist, IQ mixing, aliasing, quantisation, what the B210 delivers
     - proposal
   * - ``D_noise_sensitivity.ipynb``
     - radiometer equation, SEFD, noise of a visibility, S/N against integration time, bandwidth and baseline length
     - proposal
   * - ``E_coordinates.ipynb``
     - coordinate frames: equatorial, horizontal, WGS84/ECEF, local east-north-up, uvw
     - proposal

*available*: in the repository; *planned*: decided, not yet written;
*proposal*: suggestion for discussion. The former Appendices A (Wiener
filter vs. CLEAN) and B (the Sun as a first target) become the numbered
notebooks 07 and 08. The ``binder/`` folder holds the environment for
Binder / JupyterHub.
