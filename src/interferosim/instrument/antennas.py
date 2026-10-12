# src/interferosim/instrument/antennas.py -- antenna types (interferosim.instrument.antennas)
"""Antenna types of Interferosim.

All antennas derive from the abstract class :class:`Antenna`. Every type
provides its gain, effective area, half-power beam width (HPBW) and the
normalised power pattern (primary beam) as functions of the frequency.

Origin
------
* :class:`ParabolicDish`: ``aas6sim`` 1.3.1 ``model.Instrument``
  (``a_eff``, ``hpbw_deg``, ``u_max``, ``u_cut_eff``) and
  ``model.otf_circular``; values of notebook A (3 m, eta 0.60, k 1.20).
* :class:`Horn`: notebook B (Appendix B) and ``aas6sim.simulate.solar_transit``
  -- 15 dBi, A_eff = G lambda^2 / (4 pi) = 0.112 m^2 at 21 cm,
  HPBW = sqrt(4 pi / G / 1.133).
* :class:`Yagi`: new; the beam is derived from the gain in the same way
  as for the horn (Gaussian main lobe).

Relations (chapter "Radio Astronomical Fundamentals")
-----------------------------------------------------
* A_eff = G lambda^2 / (4 pi)
* Gaussian main lobe: Omega_A = 1.133 HPBW^2 = 4 pi / G
  =>  HPBW = sqrt(4 pi / (1.133 G))   [rad]

Version 0.3.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

import math
import tomllib
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, ClassVar, Dict, Optional, Type, Union

import numpy as np
from scipy.special import j1

from ..constants import C_LIGHT, HI_21CM

__all__ = [
    "Antenna",
    "Horn",
    "ParabolicDish",
    "Yagi",
    "ANTENNA_TYPES",
    "antenna_from_dict",
    "otf_circular",
    "hpbw_from_gain_deg",
]

# Omega_A = GAUSS_BEAM_SOLID * HPBW^2 for a Gaussian main lobe
GAUSS_BEAM_SOLID = math.pi / (4.0 * math.log(2.0))      # = 1.1331


def hpbw_from_gain_deg(gain_dbi: float) -> float:
    """HPBW in degrees of a Gaussian main lobe with the given gain
    (notebook B: HPBW = sqrt(4 pi / G / 1.133))."""
    g = 10.0 ** (gain_dbi / 10.0)
    return math.degrees(math.sqrt(4.0 * math.pi / (GAUSS_BEAM_SOLID * g)))


def otf_circular(u, u_cut: float) -> np.ndarray:
    """Transfer function of a filled circular aperture (aas6sim
    ``model.otf_circular``): autocorrelation of a uniformly illuminated
    disc, 1 at u = 0 and 0 for |u| >= u_cut."""
    q = np.clip(np.abs(np.asarray(u, dtype=float)) / u_cut, 0.0, 1.0)
    return (2.0 / np.pi) * (np.arccos(q) - q * np.sqrt(np.maximum(1.0 - q**2, 0.0)))


def _gauss_beam(theta_deg, hpbw_deg: float) -> np.ndarray:
    t = np.asarray(theta_deg, dtype=float)
    return np.exp(-4.0 * math.log(2.0) * (t / hpbw_deg) ** 2)


# ======================================================================
@dataclass
class Antenna(ABC):
    """Abstract base class of all antenna types."""

    name: str

    type: ClassVar[str] = ""

    # ------------------------------------------------------------------
    @abstractmethod
    def gain_dbi(self, f_hz: float = HI_21CM) -> float:
        """Gain over an isotropic radiator in dBi."""

    def gain(self, f_hz: float = HI_21CM) -> float:
        """Linear gain."""
        return 10.0 ** (self.gain_dbi(f_hz) / 10.0)

    def effective_area(self, f_hz: float = HI_21CM) -> float:
        """Effective area in m^2, A_eff = G lambda^2 / (4 pi)."""
        lam = C_LIGHT / f_hz
        return self.gain(f_hz) * lam**2 / (4.0 * math.pi)

    def hpbw_deg(self, f_hz: float = HI_21CM) -> float:
        """Half-power beam width in degrees (default: from the gain)."""
        return hpbw_from_gain_deg(self.gain_dbi(f_hz))

    def beam(self, theta_deg, f_hz: float = HI_21CM) -> np.ndarray:
        """Normalised power pattern (primary beam), 1 on the axis.

        ``theta_deg`` is the angle from the pointing direction. Default:
        Gaussian main lobe with :meth:`hpbw_deg` (no side lobes).
        """
        return _gauss_beam(theta_deg, self.hpbw_deg(f_hz))

    def jy_to_k(self, f_hz: float = HI_21CM) -> float:
        """Antenna temperature in K per Jy of an unresolved source
        (one polarisation), A_eff / (2 k)."""
        from ..constants import JANSKY, K_BOLTZMANN
        return self.effective_area(f_hz) * JANSKY / (2.0 * K_BOLTZMANN)

    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        """Dictionary in the format of ``configs/antennas/*.toml``."""
        d: Dict[str, Any] = {"name": self.name, "type": self.type}
        for k, v in asdict(self).items():
            if k not in d and v is not None:
                d[k] = v
        return d

    def summary(self, f_hz: float = HI_21CM) -> Dict[str, float]:
        """Key figures at one frequency (for tables)."""
        return {"gain_dbi": self.gain_dbi(f_hz),
                "a_eff_m2": self.effective_area(f_hz),
                "hpbw_deg": self.hpbw_deg(f_hz)}

    @classmethod
    def from_toml(cls, path: Union[str, Path]) -> "Antenna":
        """Read an antenna from ``configs/antennas/<name>.toml``; the key
        ``type`` selects the class."""
        with open(path, "rb") as fh:
            return antenna_from_dict(tomllib.load(fh))


# ======================================================================
@dataclass
class Horn(Antenna):
    """Horn antenna defined by its gain (kind "horn").

    Without ``design_freq_hz`` the gain is the same at all frequencies
    (as in notebook B). With ``design_freq_hz`` the horn is treated as a
    fixed aperture: G grows with f^2 and the HPBW shrinks with 1/f.
    """

    gain_dbi_nominal: float = 15.0
    design_freq_hz: Optional[float] = None

    type: ClassVar[str] = "horn"

    def gain_dbi(self, f_hz: float = HI_21CM) -> float:
        if self.design_freq_hz:
            return self.gain_dbi_nominal + 20.0 * math.log10(f_hz / self.design_freq_hz)
        return self.gain_dbi_nominal


# ======================================================================
@dataclass
class Yagi(Antenna):
    """Yagi-Uda antenna (kind "yagi"), linear polarisation.

    The gain is valid near ``design_freq_hz`` (a Yagi is narrow-band);
    the main lobe is modelled as a rotationally symmetric Gaussian derived
    from the gain (E and H plane are not distinguished).
    """

    n_elements: int = 0
    boom_length_m: float = float("nan")
    gain_dbi_nominal: float = float("nan")
    design_freq_hz: float = HI_21CM
    polarisation: str = "linear"

    type: ClassVar[str] = "yagi"

    def gain_dbi(self, f_hz: float = HI_21CM) -> float:
        return self.gain_dbi_nominal

    @property
    def boom_length_wavelengths(self) -> float:
        """Boom length in wavelengths at the design frequency."""
        return self.boom_length_m / (C_LIGHT / self.design_freq_hz)


# ======================================================================
@dataclass
class ParabolicDish(Antenna):
    """Parabolic dish (kind "dish"), values of aas6sim / notebook A.

    * A_eff = eta pi D^2 / 4
    * HPBW = beam_factor lambda / D (tapered illumination, beam_factor 1.2)
    * primary beam: Airy pattern of the equivalent uniformly illuminated
      aperture d_eff = 1.02 D / beam_factor, so that it has the same HPBW
    * transfer function :meth:`otf` with the cut-off u_cut_eff = d_eff/lambda
      (aas6sim ``u_cut_eff``; using u_max = D/lambda would make the beam too
      narrow, see the docstring there)
    """

    diameter_m: float = 3.0
    aperture_efficiency: float = 0.60
    beam_factor: float = 1.20

    type: ClassVar[str] = "dish"

    def effective_area(self, f_hz: float = HI_21CM) -> float:
        return self.aperture_efficiency * math.pi * self.diameter_m**2 / 4.0

    def gain_dbi(self, f_hz: float = HI_21CM) -> float:
        lam = C_LIGHT / f_hz
        g = 4.0 * math.pi * self.effective_area(f_hz) / lam**2
        return 10.0 * math.log10(g)

    def hpbw_deg(self, f_hz: float = HI_21CM) -> float:
        return math.degrees(self.beam_factor * (C_LIGHT / f_hz) / self.diameter_m)

    @property
    def d_eff_m(self) -> float:
        """Diameter of the uniformly illuminated aperture with the same HPBW."""
        return 1.02 * self.diameter_m / self.beam_factor

    def u_max(self, f_hz: float = HI_21CM) -> float:
        """Cut-off spatial frequency of the geometric aperture, D / lambda."""
        return self.diameter_m / (C_LIGHT / f_hz)

    def u_cut_eff(self, f_hz: float = HI_21CM) -> float:
        """Cut-off spatial frequency of the actual beam, d_eff / lambda."""
        return self.d_eff_m / (C_LIGHT / f_hz)

    def otf(self, u, f_hz: float = HI_21CM) -> np.ndarray:
        """Transfer function of the single dish for spatial frequency u [lambda]."""
        return otf_circular(u, self.u_cut_eff(f_hz))

    def beam(self, theta_deg, f_hz: float = HI_21CM) -> np.ndarray:
        x = math.pi * self.d_eff_m / (C_LIGHT / f_hz) * np.sin(np.radians(theta_deg))
        x = np.asarray(x, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            a = np.where(np.abs(x) < 1e-9, 1.0, 2.0 * j1(x) / np.where(x == 0, 1.0, x))
        return a**2


# ======================================================================
ANTENNA_TYPES: Dict[str, Type[Antenna]] = {
    Horn.type: Horn,
    Yagi.type: Yagi,
    ParabolicDish.type: ParabolicDish,
}

# keys of configs/antennas/*.toml that differ from the attribute names
_TOML_ALIASES = {"gain_dbi": "gain_dbi_nominal"}


def antenna_from_dict(d: Dict[str, Any]) -> Antenna:
    """Create an antenna from a configuration dictionary (``type`` key)."""
    kind = str(d.get("type", "")).lower()
    if kind not in ANTENNA_TYPES:
        raise ValueError(f"unknown antenna type {kind!r} (known: {', '.join(ANTENNA_TYPES)})")
    cls = ANTENNA_TYPES[kind]
    allowed = {f.name for f in fields(cls)}
    kwargs = {}
    for k, v in d.items():
        k = _TOML_ALIASES.get(k, k)
        if k in allowed:
            kwargs[k] = v
    kwargs.setdefault("name", kind)
    return cls(**kwargs)
