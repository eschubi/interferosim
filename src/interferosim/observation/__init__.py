# src/interferosim/observation/__init__.py -- sub-package interferosim.observation
"""Observation: hour angle and elevation (sidereal and solar day), uv-coverage, visibilities with noise, adding/multiplying correlator, fringe rate and delay.

Version 0.1.0: empty skeleton, filled in the migration steps.
Version 0.4.0 (2026-10-31): module geometry (GMST, source direction in ECEF,
geometric delay and rate, uvw, elevation) for the signal simulation.
Version 0.5.0 (2026-10-31): Observation, UVCoverage, VisibilitySimulator,
TransitSimulator (module obs); position of the Sun (module sun).
Version 0.8.1 (2026-10-31): the Sun as a target (module solar): sun_flux_jy,
solar_transit_mjd, SolarTransitSimulator, disc_amplitude, fit_disc_diameter.
"""

from .geometry import (SIDEREAL_DEG_PER_S, delay_rate, direction, elevation_deg,
                       geometric_delay, gha_deg, gmst_deg, sky_basis, uvw)
from .obs import Observation, TransitSimulator, UVCoverage, VisibilitySimulator
from .sun import mjd_from_iso, sun_radec_deg
from .solar import (SolarTransitSimulator, disc_amplitude, fit_disc_diameter, solar_transit_mjd,
                    sun_flux_jy)

__all__: list[str] = ["gmst_deg", "gha_deg", "sky_basis", "direction", "geometric_delay",
                      "delay_rate", "uvw", "elevation_deg", "SIDEREAL_DEG_PER_S",
                      "Observation", "UVCoverage", "VisibilitySimulator", "TransitSimulator",
                      "mjd_from_iso", "sun_radec_deg",
                      "sun_flux_jy", "solar_transit_mjd", "SolarTransitSimulator",
                      "disc_amplitude", "fit_disc_diameter"]
