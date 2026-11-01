# src/interferosim/instrument/__init__.py -- sub-package interferosim.instrument
"""Instrument: antenna types (base class + configured types), receiver model incl. Ettus USRP B210 profile, array geometry.

Version 0.1.0: empty skeleton, filled in the migration steps.
Version 0.3.0 (2026-10-12): antennas (Horn, Yagi, ParabolicDish), Receiver,
Site / Station / Array implemented.
Version 0.9.3 (2026-11-01): module sensitivity (radiometer equation, visibility
noise, integration time, image noise, coherence; theory notebook D).
"""

from .antennas import (ANTENNA_TYPES, Antenna, Horn, ParabolicDish, Yagi,
                       antenna_from_dict, hpbw_from_gain_deg, otf_circular)
from .array import (Array, Site, Station, ecef_to_enu_matrix, ecef_to_geodetic,
                    geodetic_to_ecef)
from .receiver import Receiver
from .sensitivity import (antenna_temperature_k, coherence_factor, coherence_time_s, image_noise_jy,
                          incoherent_snr, integration_time_s, phase_rms_white_fm_rad,
                          radiometer_sigma_k, random_walk_coherence, visibility_sigma_jy)

__all__: list[str] = [
    "Antenna", "Horn", "Yagi", "ParabolicDish", "ANTENNA_TYPES", "antenna_from_dict",
    "hpbw_from_gain_deg", "otf_circular",
    "Receiver",
    "Site", "Station", "Array", "geodetic_to_ecef", "ecef_to_geodetic", "ecef_to_enu_matrix",
    "antenna_temperature_k", "radiometer_sigma_k", "visibility_sigma_jy", "integration_time_s",
    "image_noise_jy", "coherence_factor", "phase_rms_white_fm_rad", "coherence_time_s",
    "random_walk_coherence", "incoherent_snr",
]
