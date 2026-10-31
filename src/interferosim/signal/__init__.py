# src/interferosim/signal/__init__.py -- sub-package interferosim.signal
"""Signal: simulated SDR receiver signals (IQ samples) of all stations,
instrument errors, quantisation, simple correlation.

Version 0.4.0 (2026-10-31): IQSimulator, IQBlock, StationErrors, quantize;
correlate (fractional_delay, correct_geometry, cross_correlate,
measure_visibilities).
"""

from .correlate import correct_geometry, cross_correlate, fractional_delay, measure_visibilities
from .iq import IQBlock, IQSimulator, StationErrors, quantize

__all__ = ["IQSimulator", "IQBlock", "StationErrors", "quantize",
           "fractional_delay", "correct_geometry", "cross_correlate", "measure_visibilities"]
