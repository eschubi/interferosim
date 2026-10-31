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
"""

from .correlate import (Accumulator, FringeSolution, apply_fringe_solution, correct_geometry,
                        cross_correlate, cross_spectrum, delay_rate_search, fractional_delay,
                        fringe_stop, integer_align, measure_visibilities, refine_fringe_solution,
                        stokes_from_products)
from .fx import FXCorrelator
from .iq import IQBlock, IQSimulator, StationErrors, quantize

__all__ = ["IQSimulator", "IQBlock", "StationErrors", "quantize",
           "fractional_delay", "correct_geometry", "cross_correlate", "measure_visibilities",
           "integer_align", "fringe_stop", "Accumulator", "cross_spectrum", "FringeSolution",
           "delay_rate_search", "apply_fringe_solution", "stokes_from_products", "FXCorrelator",
           "refine_fringe_solution"]
