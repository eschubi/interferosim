# src/interferosim/instrument/sensitivity.py -- radiometer equation, visibility noise, coherence (interferosim.instrument.sensitivity)
"""Sensitivity of a single dish and of an interferometer baseline (theory
notebook D).

* :func:`antenna_temperature_k` -- T_A = S A_eff / (2 k).
* :func:`radiometer_sigma_k` -- T_sys / sqrt(B t) (one polarisation, total
  power).
* :func:`visibility_sigma_jy` -- noise of the real (imaginary) part of one
  visibility, sqrt((SEFD_1 + S_1)(SEFD_2 + S_2) / (2 B t)) / eta_Q; S_k the
  source power received by station k (matters only for strong sources such
  as the Sun).
* :func:`integration_time_s` -- time for a given S/N of a visibility of
  amplitude ``|V|`` (vector average, the phase known or stopped).
* :func:`image_noise_jy` -- rms of a naturally weighted image from N_bl
  baselines of equal SEFD: sigma_vis / sqrt(N_bl) for the same t, B.
* :func:`coherence_factor` -- exp(-sigma_phi^2 / 2): amplitude kept by a
  coherent average over a phase with Gaussian scatter sigma_phi.
* :func:`phase_rms_white_fm_rad`, :func:`coherence_time_s` -- phase random
  walk of two independent clocks with white frequency noise,
  sigma_y(tau) = sigma_1 / sqrt(tau) each: sigma_phi(T) = 2 pi f sigma_1
  sqrt(2 T) for the difference of the two clocks; coherence time for a
  chosen phase rms.
* :func:`random_walk_coherence` -- expected amplitude of a coherent average
  over [0, T] when the phase is a random walk with sigma_phi(T) over the
  interval and the phase is referred to the middle of the interval (as a
  phase solution per segment does): 4 (1 - exp(-sigma^2/4)) / sigma^2; larger
  than exp(-sigma^2/2), because most of the interval has drifted less.
* :func:`incoherent_snr` -- detection by averaging ``|V|^2`` over n segments of
  S/N s each: S/N = s^2 sqrt(n) / (2 sqrt(1 + s^2)) (s^2 sqrt(n)/2 for weak,
  s sqrt(n)/2 for strong segments; always below the coherent s sqrt(n)).

Version 0.9.3 (2026-11-01): first implementation.
"""

from __future__ import annotations

import math

import numpy as np

from ..constants import JANSKY, K_BOLTZMANN

__all__ = ["antenna_temperature_k", "radiometer_sigma_k", "visibility_sigma_jy",
           "integration_time_s", "image_noise_jy", "coherence_factor", "phase_rms_white_fm_rad",
           "coherence_time_s", "random_walk_coherence", "incoherent_snr"]


def antenna_temperature_k(flux_jy, a_eff_m2: float):
    """Antenna temperature [K] of a source of flux density ``flux_jy`` (one
    polarisation of an unpolarised source: S/2 per polarisation is the usual
    convention; here the full S as in the rest of the package)."""
    return np.asarray(flux_jy, dtype=float) * JANSKY * a_eff_m2 / (2.0 * K_BOLTZMANN)


def radiometer_sigma_k(t_sys_k: float, bandwidth_hz: float, t_s):
    """Radiometer equation: rms of the measured system temperature."""
    return t_sys_k / np.sqrt(bandwidth_hz * np.asarray(t_s, dtype=float))


def visibility_sigma_jy(sefd_1: float, sefd_2: float, bandwidth_hz: float, t_s,
                        s_1: float = 0.0, s_2: float = 0.0, eta_q: float = 1.0):
    """Noise of the real (and of the imaginary) part of one visibility [Jy]."""
    return np.sqrt((sefd_1 + s_1) * (sefd_2 + s_2) / (2.0 * bandwidth_hz * np.asarray(t_s, dtype=float))) / eta_q


def integration_time_s(v_jy, snr: float, sefd_1: float, sefd_2: float, bandwidth_hz: float,
                       eta_q: float = 1.0, s_1: float = 0.0, s_2: float = 0.0):
    """Time for ``|V|`` / sigma = ``snr`` (sigma of the real part, see
    :func:`visibility_sigma_jy`)."""
    v = np.asarray(v_jy, dtype=float)
    return snr ** 2 * (sefd_1 + s_1) * (sefd_2 + s_2) / (2.0 * bandwidth_hz * v ** 2 * eta_q ** 2)


def image_noise_jy(sefd: float, bandwidth_hz: float, t_s: float, n_baselines: int,
                   eta_q: float = 1.0) -> float:
    """rms of a naturally weighted dirty image (point-source sensitivity) for
    ``n_baselines`` baselines of equal SEFD observing for ``t_s``."""
    return float(visibility_sigma_jy(sefd, sefd, bandwidth_hz, t_s, eta_q=eta_q)) / math.sqrt(n_baselines)


def coherence_factor(sigma_phi_rad):
    """Amplitude kept by a coherent average over a Gaussian phase scatter."""
    return np.exp(-0.5 * np.asarray(sigma_phi_rad, dtype=float) ** 2)


def phase_rms_white_fm_rad(freq_hz: float, sigma_y1: float, t_s, n_clocks: int = 2):
    """Phase rms after ``t_s`` between ``n_clocks`` independent clocks (2:
    two stations) with white FM noise sigma_y(tau) = sigma_y1 / sqrt(tau):
    sigma_phi = 2 pi f sigma_y1 sqrt(n t)."""
    return 2.0 * math.pi * freq_hz * sigma_y1 * np.sqrt(n_clocks * np.asarray(t_s, dtype=float))


def coherence_time_s(freq_hz: float, sigma_y1: float, phase_rms_rad: float = 1.0,
                     n_clocks: int = 2) -> float:
    """Time after which the phase rms of :func:`phase_rms_white_fm_rad`
    reaches ``phase_rms_rad`` (1 rad: coherence factor 0.61)."""
    return (phase_rms_rad / (2.0 * math.pi * freq_hz * sigma_y1)) ** 2 / n_clocks


def random_walk_coherence(sigma_phi_rad):
    """Expected value of (1/T) int_0^T exp(i (phi(t) - phi(T/2))) dt for a phase
    random walk whose variance grows to sigma_phi^2 over the interval T:
    4 (1 - exp(-s^2/4)) / s^2. (The magnitude of a single average with the
    best phase is slightly larger on average.)"""
    s2 = np.asarray(sigma_phi_rad, dtype=float) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(s2 < 1e-12, 1.0, 4.0 * (1.0 - np.exp(-0.25 * s2)) / np.where(s2 == 0, 1.0, s2))


def incoherent_snr(snr_segment, n_segments):
    """Detection S/N when ``|V|^2`` is averaged over ``n_segments`` segments of
    S/N ``snr_segment`` each (``|V|`` / sigma of the real part): signal ``|V|^2``,
    scatter of ``|V + noise|^2`` = 2 sigma sqrt(``|V|^2`` + sigma^2), hence
    s^2 sqrt(n) / (2 sqrt(1 + s^2))."""
    s = np.asarray(snr_segment, dtype=float)
    return s ** 2 * np.sqrt(np.asarray(n_segments, dtype=float)) / (2.0 * np.sqrt(1.0 + s ** 2))
