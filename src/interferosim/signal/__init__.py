# src/interferosim/signal/__init__.py -- sub-package interferosim.signal
"""Signal: simulated SDR receiver signals (IQ samples) of all stations,
instrument errors, quantisation, simple correlation.

Version 0.4.0 (2026-10-31): IQSimulator, IQBlock, StationErrors, quantize;
correlate (fractional_delay, correct_geometry, cross_correlate,
measure_visibilities).
Version 0.6.0 (2026-10-31): FX correlator (FXCorrelator) and the AAS6
correlation kernels (integer_align, fringe_stop, Accumulator, cross_spectrum,
FringeSolution, delay_rate_search, apply_fringe_solution, stokes_from_products),
refine_fringe_solution.
Version 0.9.1 (2026-11-01): lag_correlation, xf_cross_spectrum (XF correlator,
theory notebook B).
Version 0.9.2 (2026-11-01): module baseband (iq_imbalance, image_rejection_db,
quantizer_levels, quantization_efficiency, van_vleck_1bit; theory notebook C).
"""

from .correlate import (Accumulator, FringeSolution, apply_fringe_solution, correct_geometry,
                        cross_correlate, cross_spectrum, delay_rate_search, fractional_delay,
                        fringe_stop, integer_align, lag_correlation, measure_visibilities,
                        refine_fringe_solution, stokes_from_products, xf_cross_spectrum)
from .baseband import (image_rejection_db, iq_imbalance, quantization_efficiency,
                       quantizer_levels, van_vleck_1bit)
from .fx import FXCorrelator
from .iq import IQBlock, IQSimulator, StationErrors, quantize

__all__ = ["IQSimulator", "IQBlock", "StationErrors", "quantize",
           "fractional_delay", "correct_geometry", "cross_correlate", "measure_visibilities",
           "integer_align", "fringe_stop", "Accumulator", "cross_spectrum", "FringeSolution",
           "delay_rate_search", "apply_fringe_solution", "stokes_from_products", "FXCorrelator",
           "refine_fringe_solution", "lag_correlation", "xf_cross_spectrum",
           "iq_imbalance", "image_rejection_db", "quantizer_levels", "quantization_efficiency",
           "van_vleck_1bit"]
