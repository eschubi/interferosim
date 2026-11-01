# src/interferosim/signal/baseband.py -- complex baseband: IQ imbalance, quantisation efficiency (interferosim.signal.baseband)
"""Complex baseband and the imperfections of the receive path (theory
notebook C).

* :func:`iq_imbalance` -- apply a gain and phase error between the I and Q
  branches of a complex stream (mixer and ADC of the AD9361 in the B210);
  a tone at +f gets an image at -f.
* :func:`image_rejection_db` -- image rejection ratio
  IRR = (1 + 2 g cos(phi) + g^2) / (1 - 2 g cos(phi) + g^2), g the amplitude
  ratio of the branches, phi the phase error.
* :func:`quantizer_levels` -- thresholds and output levels of the mid-rise
  quantiser of :func:`~interferosim.signal.quantize` (same step rule).
* :func:`quantization_efficiency` -- correlation efficiency eta_Q of a
  quantiser for weakly correlated Gaussian noise (unit sigma per real
  component):  eta_Q = [sum_k y_k (phi(t_k) - phi(t_k+1))]^2 /
  sum_k y_k^2 (Phi(t_k+1) - Phi(t_k)),  phi, Phi the normal density and
  distribution; 2/pi for 1 bit, 0.881 for 2 bits at the optimum step.
* :func:`van_vleck_1bit` -- the arcsine law rho = sin(pi rho_q / 2) that
  restores the true correlation from 1-bit data.

Version 0.9.2 (2026-11-01): first implementation.
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from scipy.stats import norm

__all__ = ["iq_imbalance", "image_rejection_db", "quantizer_levels", "quantization_efficiency",
           "van_vleck_1bit", "OPTIMUM_STEP"]

# optimum step (in sigma) of the mid-rise quantiser for Gaussian noise, as in quantize()
OPTIMUM_STEP = {2: 0.996, 3: 0.586, 4: 0.335}


def iq_imbalance(x: np.ndarray, gain_db: float = 0.0, phase_deg: float = 0.0) -> np.ndarray:
    """Complex stream with the Q branch scaled by 10^(gain_db/20) and rotated
    by ``phase_deg`` against I: I' = I, Q' = g (Q cos(phi) + I sin(phi))."""
    g = 10.0 ** (gain_db / 20.0)
    p = math.radians(phase_deg)
    i, q = x.real, x.imag
    return i + 1j * g * (q * math.cos(p) + i * math.sin(p))


def image_rejection_db(gain_db: float, phase_deg: float) -> float:
    """Image rejection ratio [dB] of an IQ mixer with the given amplitude and
    phase imbalance (infinite for a perfect mixer)."""
    g = 10.0 ** (gain_db / 20.0)
    p = math.radians(phase_deg)
    num = 1.0 + 2.0 * g * math.cos(p) + g * g
    den = 1.0 - 2.0 * g * math.cos(p) + g * g
    return float("inf") if den <= 0 else 10.0 * math.log10(num / den)


def quantizer_levels(bits: int, step_sigma: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """(thresholds, outputs) of the mid-rise quantiser in units of sigma:
    outputs (k + 1/2) step for k = -2^(b-1) .. 2^(b-1) - 1, thresholds between
    them (the outer ones at -inf / +inf). 1 bit: outputs -1, +1, threshold 0.
    Default step as in :func:`~interferosim.signal.quantize`."""
    if bits < 1:
        raise ValueError("bits must be >= 1")
    if bits == 1:
        return np.array([-np.inf, 0.0, np.inf]), np.array([-1.0, 1.0])
    full = 2 ** (bits - 1)
    step = step_sigma if step_sigma is not None else OPTIMUM_STEP.get(bits, 1.0 / (0.25 * full))
    k = np.arange(-full, full)
    outputs = (k + 0.5) * step
    inner = (np.arange(-full + 1, full)) * step
    return np.concatenate([[-np.inf], inner, [np.inf]]), outputs


def quantization_efficiency(bits: int, step_sigma: Optional[float] = None) -> float:
    """Correlation efficiency eta_Q (weak signals) of the quantiser of
    :func:`quantizer_levels`; 1 for an infinite number of levels."""
    t, y = quantizer_levels(bits, step_sigma)
    pdf = norm.pdf(t)                      # 0 at +-inf
    cdf = norm.cdf(t)
    num = np.sum(y * (pdf[:-1] - pdf[1:])) ** 2
    den = np.sum(y ** 2 * (cdf[1:] - cdf[:-1]))
    return float(num / den)


def van_vleck_1bit(rho_q) -> np.ndarray:
    """True correlation coefficient from the correlation of 1-bit (sign)
    data: rho = sin(pi rho_q / 2) (Van Vleck, arcsine law)."""
    return np.sin(0.5 * np.pi * np.asarray(rho_q, dtype=float))
