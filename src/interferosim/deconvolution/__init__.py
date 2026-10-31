# src/interferosim/deconvolution/__init__.py -- sub-package interferosim.deconvolution
"""Deconvolution: Hoegbom-CLEAN, Wiener filter, clean beam and restoration.

Version 0.1.0: empty skeleton, filled in the migration steps.
Version 0.8.0 (2026-10-31): Deconvolver, CleanBeam, HogbomClean, DeconvResult, WienerFilter,
reference_image, rms_error (module core).
"""

from .core import CleanBeam, DeconvResult, Deconvolver, HogbomClean, WienerFilter, reference_image, rms_error

__all__: list[str] = ["Deconvolver", "CleanBeam", "HogbomClean", "DeconvResult", "WienerFilter",
                      "reference_image", "rms_error"]
