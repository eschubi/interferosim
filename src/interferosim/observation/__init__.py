# src/interferosim/observation/__init__.py -- sub-package interferosim.observation
"""Observation: hour angle and elevation (sidereal and solar day), uv-coverage, visibilities with noise, adding/multiplying correlator, fringe rate and delay.

Version 0.1.0: empty skeleton, filled in the migration steps.
Version 0.4.0 (2026-10-31): module geometry (GMST, source direction in ECEF,
geometric delay and rate, uvw, elevation) for the signal simulation.
"""

from .geometry import (SIDEREAL_DEG_PER_S, delay_rate, direction, elevation_deg,
                       geometric_delay, gha_deg, gmst_deg, sky_basis, uvw)

__all__: list[str] = ["gmst_deg", "gha_deg", "sky_basis", "direction", "geometric_delay",
                      "delay_rate", "uvw", "elevation_deg", "SIDEREAL_DEG_PER_S"]
