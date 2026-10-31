.. docs/sphinx/interferometer.rst -- the radio interferometer in brief

The radio interferometer
========================

A single radio telescope of diameter :math:`D` resolves angles of about
:math:`\lambda / D` -- at the 21 cm line a 3 m dish sees a blur of several
degrees. An interferometer combines the signals of two or more antennas a
distance :math:`B` apart (the *baseline*) and reaches a resolution of about
:math:`\lambda / B`. Interferosim simulates this chain from the sky to the
image.

Two-element interferometer
--------------------------

A plane wave from the direction :math:`\hat{s}` reaches the second antenna
later than the first by the geometric delay

.. math::

   \tau_g = \frac{\vec{B} \cdot \hat{s}}{c}.

The correlator multiplies and averages the two voltages. For a point source
of flux density :math:`S` the output oscillates with the source position --
the *fringe* pattern:

.. math::

   r(\tau_g) \propto S \, \cos(2 \pi \nu \tau_g).

As the earth rotates, :math:`\hat{s}` changes and the fringes run through the
correlator output; their frequency grows with the baseline length.

Visibility and uv-plane
-----------------------

The baseline measured in wavelengths and projected onto the plane
perpendicular to the field centre gives the coordinates :math:`(u, v)`.
The complex correlator output is the *visibility*. By the van
Cittert--Zernike theorem it is the Fourier transform of the sky brightness
:math:`I(l, m)` (direction cosines relative to the field centre), weighted
with the primary beam :math:`A(l, m)` of the antennas:

.. math::

   V(u, v) = \iint A(l, m)\, I(l, m)\,
             e^{-2 \pi i (u l + v m)} \, \mathrm{d}l \, \mathrm{d}m.

Each baseline samples one point :math:`(u, v)` per instant (and its mirror
point :math:`(-u, -v)`); earth rotation draws tracks in the uv-plane. Short
baselines see extended structure, long baselines fine detail -- a source much
larger than :math:`\lambda / B` is "resolved out". The sky models of
Interferosim provide :math:`V(u, v)` analytically
(:meth:`interferosim.sky.SkyModel.visibility`).

Sensitivity
-----------

The noise depends on the system temperature :math:`T_\mathrm{sys}`, the
effective area :math:`A_\mathrm{eff}` of the antennas, the bandwidth
:math:`B_\nu` and the integration time :math:`\tau`. With the system
equivalent flux density

.. math::

   \mathrm{SEFD} = \frac{2 k T_\mathrm{sys}}{A_\mathrm{eff}}

the noise of one baseline :math:`i`--:math:`j` (real and imaginary part
each) is

.. math::

   \sigma = \frac{\sqrt{\mathrm{SEFD}_i \, \mathrm{SEFD}_j}}{\sqrt{2 B_\nu \tau}}.

:class:`interferosim.instrument.Receiver` computes both
(``sefd_jy``, ``baseline_sigma_jy``).

From visibilities to an image
-----------------------------

The sampled visibilities are gridded, weighted and inverse Fourier
transformed into the *dirty image* -- the true sky convolved with the
*dirty beam* (the point spread function of the uv-coverage).
Deconvolution (Högbom-CLEAN, Wiener filter) removes the sidelobes of the
dirty beam and restores the image.

Processing chain of Interferosim
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Stage
     - Module
   * - Sky models (point, extended, HI 21 cm, OH masers, solar disc)
     - :mod:`interferosim.sky`
   * - Antennas (selectable types), receiver (Ettus USRP B210), array
     - :mod:`interferosim.instrument`
   * - Simulated SDR receiver signals (IQ), instrument errors
     - :mod:`interferosim.signal`
   * - Geometry, uv-coverage, visibilities, correlator
     - :mod:`interferosim.observation` (geometry, uv-coverage, visibility and transit simulation; correlator in preparation), :mod:`interferosim.visdata`
   * - Gridding, weighting, dirty image and beam
     - :mod:`interferosim.imaging`
   * - Högbom-CLEAN, Wiener filter, restoration
     - :mod:`interferosim.deconvolution` (in preparation)
   * - FITS input/output
     - :mod:`interferosim.io` (in preparation)

The full derivations are in the LaTeX documentation (``docs/interferosim.tex``,
part "Theory").
