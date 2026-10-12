# src/interferosim/constants.py -- physical constants and spectral lines
"""Physical constants (SI) and rest frequencies of the simulated lines.

Values as used in the teaching notebooks (Appendix A/B) and the project
description.
"""

C_LIGHT = 299_792_458.0          # m/s
K_BOLTZMANN = 1.380649e-23       # J/K
JANSKY = 1e-26                   # W m^-2 Hz^-1
SFU = 1e4 * JANSKY               # solar flux unit

OMEGA_SIDEREAL = 7.292115e-5     # rad/s, Earth rotation (sidereal day)
OMEGA_SOLAR = 2.0 * 3.141592653589793 / 86_400.0   # rad/s, solar day

# Rest frequencies in Hz
HI_21CM = 1_420_405_751.768
OH_1612 = 1_612.231e6            # satellite line
OH_1665 = 1_665.402e6            # main line
OH_1667 = 1_667.359e6            # main line
OH_1720 = 1_720.530e6            # satellite line

LINES = {
    "HI": HI_21CM,
    "OH1612": OH_1612,
    "OH1665": OH_1665,
    "OH1667": OH_1667,
    "OH1720": OH_1720,
}


def wavelength(frequency_hz: float) -> float:
    """Wavelength in metres for a frequency in Hz."""
    return C_LIGHT / frequency_hz
