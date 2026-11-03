# src/interferosim/__init__.py -- package interferosim
"""Interferosim -- open-source simulator for radio interferometry.

Sub-packages (filled step by step from aas6sim 1.1.0 and the teaching
notebooks):

sky            sky models
instrument     antennas, receiver, array
signal         simulated SDR receiver signals (IQ), correlation check
observation    geometry, uv-coverage, visibilities, correlator
imaging        gridding, weighting, dirty image and beam
deconvolution  Hoegbom-CLEAN, Wiener filter, restoration
io             FITS input/output
visdata        visibility data class VisData (module)
plotting       reusable plot functions
widgets        ipywidgets front end
"""

__version__ = "1.0.0"

__all__ = ["__version__"]
