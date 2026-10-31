# src/interferosim/sky/spectra.py -- spectra of sky components (interferosim.sky.spectra)
"""Spectra of sky components: continuum and spectral lines.

A spectrum gives the *shape* of a component over frequency: the flux
density at the frequency f (in the LSR frame) is ``flux_jy * shape(f)``.
For a continuum source ``flux_jy`` is the flux density at the reference
frequency, for an emission line (maser) the peak flux density of the line.

* :class:`FlatSpectrum` -- shape 1 (the default of every component).
* :class:`PowerLawSpectrum` -- (f / f_ref)^alpha (synchrotron: alpha < 0).
* :class:`GaussianLineSpectrum` -- Gaussian line in velocity, optional
  continuum level relative to the peak.
* :class:`AbsorptionSpectrum` -- continuum (power law) times exp(-tau(v)),
  tau a sum of Gaussian components (HI absorption against a background
  source).

Velocities use the **radio convention**, v = c (f_0 - f) / f_0, relative to
the rest frequency f_0 of the line (:func:`radio_velocity_kms`,
:func:`frequency_from_velocity`). In the JSON sky model a component may carry
a key ``"spectrum"`` with ``{"type": ..., parameters}`` (:func:`spectrum_from_dict`).

Version 0.9.0 (2026-11-01): first implementation.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Sequence, Tuple, Type

import numpy as np

C_KMS = 299_792.458

__all__ = ["Spectrum", "FlatSpectrum", "PowerLawSpectrum", "GaussianLineSpectrum",
           "AbsorptionSpectrum", "spectrum_from_dict", "SPECTRUM_TYPES",
           "radio_velocity_kms", "frequency_from_velocity", "C_KMS"]


def radio_velocity_kms(freq_hz, rest_hz: float) -> np.ndarray:
    """Radio velocity v = c (f0 - f) / f0 in km/s."""
    return C_KMS * (1.0 - np.asarray(freq_hz, dtype=float) / rest_hz)


def frequency_from_velocity(v_kms, rest_hz: float) -> np.ndarray:
    """Frequency of the radio velocity ``v_kms``: f = f0 (1 - v/c)."""
    return rest_hz * (1.0 - np.asarray(v_kms, dtype=float) / C_KMS)


def _gauss(v, v0, fwhm):
    return np.exp(-4.0 * math.log(2.0) * ((v - v0) / fwhm) ** 2)


# ======================================================================
class Spectrum(ABC):
    """Shape of a component spectrum; flux density = flux_jy * shape(f)."""

    type: str = ""

    @abstractmethod
    def shape(self, freq_hz) -> np.ndarray:
        """Relative flux density at the (LSR-frame) frequencies ``freq_hz``."""

    def to_dict(self) -> Dict[str, Any]:
        d = {"type": self.type}
        d.update(asdict(self))
        return d


@dataclass
class FlatSpectrum(Spectrum):
    """Frequency-independent flux density."""

    type = "flat"

    def shape(self, freq_hz) -> np.ndarray:
        return np.ones_like(np.asarray(freq_hz, dtype=float))


@dataclass
class PowerLawSpectrum(Spectrum):
    """S(f) = flux_jy (f / ref_hz)^alpha."""

    alpha: float = -0.7
    ref_hz: float = 1.420405751768e9
    type = "power_law"

    def shape(self, freq_hz) -> np.ndarray:
        return (np.asarray(freq_hz, dtype=float) / self.ref_hz) ** self.alpha


@dataclass
class GaussianLineSpectrum(Spectrum):
    """Gaussian emission line at the LSR velocity ``v_kms`` with full width
    ``fwhm_kms``; ``continuum`` is a flat continuum relative to the peak."""

    rest_hz: float = 1.665402e9
    v_kms: float = 0.0
    fwhm_kms: float = 1.0
    continuum: float = 0.0
    type = "gauss_line"

    def shape(self, freq_hz) -> np.ndarray:
        v = radio_velocity_kms(freq_hz, self.rest_hz)
        return self.continuum + _gauss(v, self.v_kms, self.fwhm_kms)


@dataclass
class AbsorptionSpectrum(Spectrum):
    """Continuum (power law, ``alpha`` at ``ref_hz``) attenuated by
    exp(-tau(v)); tau(v) = sum of Gaussians given as (v_kms, fwhm_kms, tau0)."""

    rest_hz: float = 1.420405751768e9
    lines: List[Tuple[float, float, float]] = field(default_factory=list)
    alpha: float = 0.0
    ref_hz: float = 1.420405751768e9
    type = "absorption"

    def tau(self, freq_hz) -> np.ndarray:
        v = radio_velocity_kms(freq_hz, self.rest_hz)
        t = np.zeros_like(v)
        for v0, fwhm, tau0 in self.lines:
            t = t + tau0 * _gauss(v, v0, fwhm)
        return t

    def shape(self, freq_hz) -> np.ndarray:
        f = np.asarray(freq_hz, dtype=float)
        return (f / self.ref_hz) ** self.alpha * np.exp(-self.tau(f))

    def to_dict(self) -> Dict[str, Any]:
        d = super().to_dict()
        d["lines"] = [list(x) for x in self.lines]
        return d


SPECTRUM_TYPES: Dict[str, Type[Spectrum]] = {
    FlatSpectrum.type: FlatSpectrum,
    PowerLawSpectrum.type: PowerLawSpectrum,
    GaussianLineSpectrum.type: GaussianLineSpectrum,
    AbsorptionSpectrum.type: AbsorptionSpectrum,
}


def spectrum_from_dict(d: Dict[str, Any]) -> Spectrum:
    """Spectrum from its JSON form ``{"type": ..., parameters}``."""
    t = str(d.get("type", "flat"))
    if t not in SPECTRUM_TYPES:
        raise ValueError(f"unknown spectrum type {t!r} (known: {', '.join(SPECTRUM_TYPES)})")
    kw = {k: v for k, v in d.items() if k != "type"}
    if t == "absorption" and "lines" in kw:
        kw["lines"] = [tuple(float(y) for y in x) for x in kw["lines"]]
    return SPECTRUM_TYPES[t](**kw)
