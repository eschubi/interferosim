# src/interferosim/signal/correlate.py -- delay correction and simple correlation (interferosim.signal.correlate)
"""Minimal correlation of simulated (or recorded) IQ streams.

Enough to check the signal simulation against the analytic visibilities;
the full FX correlator (cross spectra, delay/rate search) follows with
notebook 05.

* :func:`fractional_delay` -- delay a time series by a (non-integer)
  number of samples with a linear phase ramp in the frequency domain
  (AAS6 ``modules/interfero/correlate.py`` 1.0.0, unchanged).
* :func:`correct_geometry` -- remove the geometric delay and the fringe
  phase of the phase centre from one station's stream (delay tracking and
  fringe stopping with the model delay).
* :func:`cross_correlate` -- average x_i x_j* over intervals.
* :func:`measure_visibilities` -- all baselines of an
  :class:`~interferosim.signal.IQSimulator` run.

Version 0.4.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

from typing import Callable, Dict, Optional

import numpy as np

__all__ = ["fractional_delay", "correct_geometry", "cross_correlate", "measure_visibilities"]


def fractional_delay(x: np.ndarray, delay_samples: float) -> np.ndarray:
    """Delay x by ``delay_samples`` (may be non-integer). Positive delay:
    the signal arrives LATER (shifted to the right). Linear phase ramp in
    the frequency domain (all-pass, circular -- the caller provides
    enough margin)."""
    x = np.asarray(x)
    n = x.shape[-1]
    k = np.fft.fftfreq(n)                       # cycles per sample
    ramp = np.exp(-2j * np.pi * k * delay_samples)
    return np.fft.ifft(np.fft.fft(x, axis=-1) * ramp, axis=-1)


def correct_geometry(x: np.ndarray, delay_s: Callable[[np.ndarray], np.ndarray],
                     fs: float, fc: float, segment: int = 16384,
                     lo_offset_hz: float = 0.0) -> np.ndarray:
    """Undo the geometric delay tau(t) of one station: advance the stream by
    tau (fractional delay -tau, per segment at its mid time) and multiply by
    exp(+2 pi i f_c tau(t)) (fringe stopping). ``delay_s(t)`` gives the
    model delay at the times t [s]; ``lo_offset_hz`` an assumed LO offset to
    remove as well."""
    x = np.asarray(x)
    out = np.empty(x.shape, np.complex128)
    n_total = x.size
    for b0 in range(0, n_total, segment):
        n = min(segment, n_total - b0)
        t = (b0 + np.arange(n)) / fs
        t_mid = t[0] + 0.5 * n / fs
        tau_mid = float(delay_s(np.array([t_mid]))[0])
        # margin from the neighbouring samples against the circular shift
        pad = int(np.ceil(abs(tau_mid) * fs)) + 8
        lo, hi = max(b0 - pad, 0), min(b0 + n + pad, n_total)
        seg = fractional_delay(x[lo:hi], -tau_mid * fs)[b0 - lo:b0 - lo + n]
        tau_t = delay_s(t)
        ph = 2 * np.pi * fc * tau_t - 2 * np.pi * lo_offset_hz * t
        out[b0:b0 + n] = seg * np.exp(1j * ph)
    return out


def cross_correlate(a: np.ndarray, b: np.ndarray, n_avg: int) -> np.ndarray:
    """Average a b* over consecutive intervals of ``n_avg`` samples."""
    n = (min(a.size, b.size) // n_avg) * n_avg
    p = a[:n] * np.conj(b[:n])
    return p.reshape(-1, n_avg).mean(axis=1)


def measure_visibilities(blocks: Dict[str, "object"], sim, t_avg_s: float,
                         lha0_deg: float = 0.0, correct: bool = True,
                         lo_offsets_hz: Optional[Dict[str, float]] = None) -> Dict[str, dict]:
    """Correlate all baselines of a simulation run.

    Returns {"A-B": {"t_s", "uvw" (n_t, 3), "vis" (n_t,), "sigma"}} with one
    entry per averaging interval of ``t_avg_s``. ``correct=False`` skips
    delay tracking and fringe stopping (raw fringes); ``lo_offsets_hz``
    removes assumed LO offsets per station (e.g. a fringe-rate solution).
    """
    st = sim.array.stations
    names = [s.name for s in st]
    fs, fc = sim.fs, sim.fc
    n_avg = max(int(round(t_avg_s * fs)), 1)
    lo_offsets_hz = lo_offsets_hz or {}
    streams = []
    for k, name in enumerate(names):
        x = blocks[name].data
        if correct and k > 0:
            x = correct_geometry(x, lambda t, k=k: sim.phase_centre_delay_s(k, t, lha0_deg),
                                 fs, fc, lo_offset_hz=lo_offsets_hz.get(name, 0.0))
        elif lo_offsets_hz.get(name):
            t = np.arange(x.size) / fs
            x = x * np.exp(-2j * np.pi * lo_offsets_hz[name] * t)
        streams.append(x)
    res = {}
    for (i, j), bl in zip(sim.array.baseline_pairs(), sim.array.baseline_names()):
        vis = cross_correlate(streams[i], streams[j], n_avg)
        t = (np.arange(vis.size) + 0.5) * n_avg / fs
        u = sim.uvw(t, lha0_deg)[sim.array.baseline_pairs().index((i, j))]
        res[bl] = {"t_s": t, "uvw": u, "vis": vis, "i": i, "j": j,
                   "sigma": sim.sigma_visibility_jy(i, j, n_avg / fs)}
    return res
