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

Correlation kernels taken over from AAS6 ``modules/interfero/correlate.py``
1.0.0 (comments translated, behaviour unchanged), since 0.6.0:

* :func:`integer_align` -- shift one stream by whole samples and crop both.
* :func:`fringe_stop` -- remove a linear fringe (phase) drift.
* :class:`Accumulator`, :func:`cross_spectrum` -- mean cross spectrum
  X_A X_B* over many FFT windows, auto spectra, coherence.
* :class:`FringeSolution`, :func:`delay_rate_search`,
  :func:`apply_fringe_solution` -- 2-D search of residual delay and residual
  fringe rate on a bright segment, and its application.
* :func:`stokes_from_products` -- Stokes visibilities from the four
  station cross products (linear or circular basis).

New in Interferosim 0.6.0: :func:`refine_fringe_solution` -- refines a
solution of :func:`delay_rate_search` by fitting the phase slope across the
band (delay) and the phase drift in time (rate); the quadratic peak
interpolation of the search has a sub-sample bias of up to ~0.3 samples.

New in Interferosim 0.9.1 (theory notebook B): :func:`lag_correlation` --
the cross-correlation function R(k) = <a[n+k] b*[n]> over the lags
-n/2 .. n/2-1, and :func:`xf_cross_spectrum` -- the cross spectrum of the XF
correlator (correlate first, then Fourier transform the lags; Wiener-Khinchin).
Its channels have the layout and the unit of ``Accumulator.visibility() /
window_power`` with a rectangular window.

All spectral averages are linear in complex amplitude or power -- never in
dB (AAS6 rule).

Version 0.4.0 (2026-10-31): first implementation.
Version 0.6.0 (2026-10-31): AAS6 correlation kernels added (see above).
Version 0.9.1 (2026-11-01): lag_correlation, xf_cross_spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

import numpy as np

__all__ = ["fractional_delay", "correct_geometry", "cross_correlate", "measure_visibilities",
           "integer_align", "fringe_stop", "Accumulator", "cross_spectrum", "FringeSolution",
           "delay_rate_search", "apply_fringe_solution", "stokes_from_products",
           "refine_fringe_solution", "lag_correlation", "xf_cross_spectrum"]


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


# ======================================================================
# AAS6 correlation kernels (modules/interfero/correlate.py 1.0.0)
# ======================================================================
def integer_align(a: np.ndarray, b: np.ndarray, shift_b: int) -> Tuple[np.ndarray, np.ndarray]:
    """Shift b by ``shift_b`` WHOLE samples against a and crop both to the
    common valid range. Positive shift_b: b is pulled EARLIER
    (b[shift_b:] against a[:-shift_b])."""
    a = np.asarray(a)
    b = np.asarray(b)
    if shift_b == 0:
        n = min(a.shape[-1], b.shape[-1])
        return a[..., :n], b[..., :n]
    if shift_b > 0:
        b2 = b[..., shift_b:]
        n = min(a.shape[-1], b2.shape[-1])
        return a[..., :n], b2[..., :n]
    s = -shift_b
    a2 = a[..., s:]
    n = min(a2.shape[-1], b.shape[-1])
    return a2[..., :n], b[..., :n]


def fringe_stop(x: np.ndarray, fringe_hz: float, fs: float, t0: float = 0.0,
                phase0: float = 0.0) -> np.ndarray:
    """Multiply x by exp(-j (2 pi fringe_hz t + phase0)), t from t0 in seconds
    (sample spacing 1/fs). Removes a linear phase (fringe) drift between the
    stations."""
    x = np.asarray(x)
    n = x.shape[-1]
    t = t0 + np.arange(n) / fs
    return x * np.exp(-1j * (2.0 * np.pi * fringe_hz * t + phase0))


@dataclass
class Accumulator:
    """Accumulates the mean cross spectrum over many FFT windows."""

    nfft: int
    window: np.ndarray
    cross: np.ndarray            # complex sum X_A conj(X_B)
    auto_a: np.ndarray           # sum |X_A|^2
    auto_b: np.ndarray           # sum |X_B|^2
    n_win: int = 0

    @classmethod
    def create(cls, nfft: int, window: str = "hann") -> "Accumulator":
        if window == "hann":
            w = np.hanning(nfft).astype(np.float64)
        elif window in ("rect", "boxcar", None):
            w = np.ones(nfft, dtype=np.float64)
        else:
            raise ValueError("unknown window: %s" % window)
        z = np.zeros(nfft, dtype=np.complex128)
        return cls(nfft=nfft, window=w, cross=z.copy(),
                   auto_a=np.zeros(nfft), auto_b=np.zeros(nfft))

    def add(self, a: np.ndarray, b: np.ndarray) -> int:
        """a, b: complex time series of equal length (already aligned and
        fringe-stopped). Split into non-overlapping nfft windows, FFT, add
        X_A conj(X_B). Returns the number of new windows."""
        a = np.asarray(a); b = np.asarray(b)
        n = min(a.shape[-1], b.shape[-1])
        nwin = n // self.nfft
        if nwin == 0:
            return 0
        aa = a[:nwin * self.nfft].reshape(nwin, self.nfft) * self.window
        bb = b[:nwin * self.nfft].reshape(nwin, self.nfft) * self.window
        FA = np.fft.fftshift(np.fft.fft(aa, axis=1), axes=1)
        FB = np.fft.fftshift(np.fft.fft(bb, axis=1), axes=1)
        self.cross += (FA * np.conj(FB)).sum(axis=0)
        self.auto_a += (np.abs(FA) ** 2).sum(axis=0)
        self.auto_b += (np.abs(FB) ** 2).sum(axis=0)
        self.n_win += nwin
        return nwin

    # -- results -------------------------------------------------------
    def visibility(self) -> np.ndarray:
        """Mean cross spectrum (complex, per frequency channel)."""
        if self.n_win == 0:
            return np.zeros(self.nfft, dtype=np.complex128)
        return self.cross / self.n_win

    def coherence(self) -> np.ndarray:
        """Normalised magnitude ``|<X_A conj(X_B)>| / sqrt(<|X_A|^2><|X_B|^2>)``
        per channel -- 0..1, robust against level differences of the stations."""
        if self.n_win == 0:
            return np.zeros(self.nfft)
        denom = np.sqrt(self.auto_a * self.auto_b) + 1e-30
        return np.abs(self.cross) / denom

    def freq_axis(self, fs: float, center_hz: float = 0.0) -> np.ndarray:
        return center_hz + np.fft.fftshift(np.fft.fftfreq(self.nfft, d=1.0 / fs))

    # -- added in Interferosim 0.6.0 -------------------------------------
    @property
    def window_power(self) -> float:
        """sum w^2 -- the cross spectrum divided by it is in the unit of the
        per-sample correlation (Jy for the simulated streams)."""
        return float(np.sum(self.window ** 2))


def cross_spectrum(a: np.ndarray, b: np.ndarray, nfft: int,
                   window: str = "hann") -> Tuple[np.ndarray, np.ndarray, int]:
    """Convenience: (visibility, coherence, n_win) for a, b."""
    acc = Accumulator.create(nfft, window)
    acc.add(a, b)
    return acc.visibility(), acc.coherence(), acc.n_win


@dataclass
class FringeSolution:
    delay_samples: float
    rate_hz: float
    peak: float                  # coherence peak 0..1
    snr: float                   # peak / median background of the search surface
    delay_grid: np.ndarray
    rate_grid: np.ndarray
    surface: np.ndarray          # |correlation| (len(rate_grid), len(delay_grid))


def delay_rate_search(a: np.ndarray, b: np.ndarray, fs: float,
                      max_delay_samples: float = 64.0,
                      delay_step: float = 0.25,
                      max_rate_hz: float = 5.0,
                      n_rate: int = 201) -> FringeSolution:
    """2-D search of residual delay and residual fringe rate between a and b.

    For every trial rate b is fringe-stopped, then the delay is found as the
    position of the envelope peak of the cross power spectrum (FFT of
    a conj(b)). Returns the best (delay_samples, rate_hz) and the whole search
    surface for diagnosis.

    Meant for a short bright segment (some 0.1-1 s); the solution then feeds
    fractional_delay()/fringe_stop() for the whole track."""
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(b, dtype=np.complex128)
    n = min(a.shape[-1], b.shape[-1])
    a = a[:n]; b = b[:n]

    rates = np.linspace(-max_rate_hz, max_rate_hz, int(n_rate))
    t = np.arange(n) / fs
    # delay axis from the FFT length; only the window of interest
    lag = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / n))     # -n/2 .. n/2-1
    sel = np.abs(lag) <= max_delay_samples
    lag_sel = lag[sel]

    surface = np.zeros((rates.size, lag_sel.size))
    A = np.fft.fft(a)
    for i, r in enumerate(rates):
        bs = b * np.exp(-2j * np.pi * r * t)
        xps = A * np.conj(np.fft.fft(bs))
        corr = np.fft.fftshift(np.fft.ifft(xps))
        surface[i] = np.abs(corr[sel])

    # finer in delay: quadratic interpolation around the maximum
    ri, di = np.unravel_index(np.argmax(surface), surface.shape)
    peak = float(surface[ri, di])
    bg = float(np.median(surface)) + 1e-30
    # sub-sample refinement along the delay
    d_hat = float(lag_sel[di])
    if 0 < di < lag_sel.size - 1:
        y0, y1, y2 = surface[ri, di - 1], surface[ri, di], surface[ri, di + 1]
        denom = (y0 - 2 * y1 + y2)
        if abs(denom) > 1e-30:
            d_hat += 0.5 * (y0 - y2) / denom * (lag_sel[1] - lag_sel[0])
    # optionally along the rate
    r_hat = float(rates[ri])
    if 0 < ri < rates.size - 1:
        y0, y1, y2 = surface[ri - 1, di], surface[ri, di], surface[ri + 1, di]
        denom = (y0 - 2 * y1 + y2)
        if abs(denom) > 1e-30:
            r_hat += 0.5 * (y0 - y2) / denom * (rates[1] - rates[0])

    # rough normalisation of the coherence peak (energy of both series)
    norm = np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2)) + 1e-30
    coh_peak = peak / norm * n            # heuristic normalisation

    return FringeSolution(delay_samples=d_hat, rate_hz=r_hat,
                          peak=min(coh_peak, 1.0), snr=peak / bg,
                          delay_grid=lag_sel / fs, rate_grid=rates,
                          surface=surface)


def apply_fringe_solution(b: np.ndarray, sol: "FringeSolution", fs: float,
                          t0: float = 0.0) -> np.ndarray:
    """Apply the solution of delay_rate_search() to b, so that b is then
    ALIGNED with a (common phase centre).

    Convention (consistent with the kernels of this module, verified by test):

    - delay_samples is the lag position of the cross-correlation peak;
      fractional_delay(b, delay_samples) shifts b exactly so that the peak
      moves to zero.
    - rate_hz is the residual fringe rate; fringe_stop(b, rate_hz, fs)
      removes it.

    Order: first the delay, then fringe stopping."""
    return fringe_stop(fractional_delay(np.asarray(b), sol.delay_samples),
                       sol.rate_hz, fs, t0=t0)


def stokes_from_products(v_xx: np.ndarray, v_yy: Optional[np.ndarray],
                         v_xy: Optional[np.ndarray] = None,
                         v_yx: Optional[np.ndarray] = None,
                         basis: str = "linear") -> dict:
    """Stokes visibilities from the station cross products.

    v_xx = <A_ch0 conj(B_ch0)>,  v_yy = <A_ch1 conj(B_ch1)>,
    v_xy = <A_ch0 conj(B_ch1)>,  v_yx = <A_ch1 conj(B_ch0)>.

    basis "linear" (X/Y dipoles):
        I = v_xx + v_yy,  Q = v_xx - v_yy,  U = v_xy + v_yx,  V = 1j (v_yx - v_xy)
    basis "circular" (R/L):
        I = v_rr + v_ll,  V = v_rr - v_ll,  Q = v_rl + v_lr,  U = 1j (v_lr - v_rl)

    With only one channel per station (v_yy/v_xy/v_yx None) only I is
    returned (I = v_xx)."""
    out = {}
    if v_yy is None:
        out["I"] = np.asarray(v_xx)
        return out
    v_xx = np.asarray(v_xx); v_yy = np.asarray(v_yy)
    if v_xy is None or v_yx is None:
        out["I"] = v_xx + v_yy
        out["Q"] = v_xx - v_yy
        return out
    v_xy = np.asarray(v_xy); v_yx = np.asarray(v_yx)
    if basis == "linear":
        out["I"] = v_xx + v_yy
        out["Q"] = v_xx - v_yy
        out["U"] = v_xy + v_yx
        out["V"] = 1j * (v_yx - v_xy)
    elif basis == "circular":
        out["I"] = v_xx + v_yy
        out["V"] = v_xx - v_yy
        out["Q"] = v_xy + v_yx
        out["U"] = 1j * (v_yx - v_xy)
    else:
        raise ValueError("basis must be 'linear' or 'circular'")
    return out


# ======================================================================
# Interferosim 0.6.0: refinement of a fringe solution
# ======================================================================
def refine_fringe_solution(a: np.ndarray, b: np.ndarray, sol: FringeSolution, fs: float,
                           nfft: int = 64, n_seg: int = 8, iterations: int = 2) -> FringeSolution:
    """Refine delay and rate of ``sol`` on the same segment (a, b).

    After applying the solution, the remaining delay eps shows as the phase
    slope 2 pi f eps of the cross spectrum a b'* across the band, the
    remaining rate as a phase drift of the channel-averaged correlation over
    ``n_seg`` sub-segments. Both are fitted (weighted by the coherence) and
    added to the solution; repeated ``iterations`` times. Returns a new
    FringeSolution (search surface of the original kept)."""
    from dataclasses import replace as _replace
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(b, dtype=np.complex128)
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    out = sol
    f = np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / fs))
    for _ in range(iterations):
        bb = apply_fringe_solution(b, out, fs)
        vis, coh, nw = cross_spectrum(a, bb, nfft)
        if nw == 0:
            break
        w = coh ** 2
        ph = np.unwrap(np.angle(vis))
        slope = np.polyfit(f, ph, 1, w=np.sqrt(w))[0]              # rad / Hz
        eps = slope / (2.0 * np.pi)                                  # s: b' delayed by eps
        new_delay = out.delay_samples - eps * fs
        out = _replace(out, delay_samples=float(new_delay))
        # rate from the phase drift over sub-segments
        bb = apply_fringe_solution(b, out, fs)
        m = n // n_seg
        if m < nfft:
            continue
        seg_ph = np.unwrap([np.angle(np.mean(a[k * m:(k + 1) * m] * np.conj(bb[k * m:(k + 1) * m])))
                            for k in range(n_seg)])
        t = (np.arange(n_seg) + 0.5) * m / fs
        drift = np.polyfit(t, seg_ph, 1)[0]                          # rad / s
        out = _replace(out, rate_hz=float(out.rate_hz - drift / (2.0 * np.pi)))
    return out


# -- added in Interferosim 0.9.1 (theory notebook B) ------------------------
def lag_correlation(a: np.ndarray, b: np.ndarray, n_lags: int) -> Tuple[np.ndarray, np.ndarray]:
    """Cross-correlation function R(k) = mean_n a[n+k] conj(b[n]) for the
    lags k = -n_lags/2 .. n_lags/2 - 1 (the same samples n for every lag).
    Returns (lags, R). If b is a delayed by d samples (b[n] = a[n-d]), the
    peak lies at k = -d; if a is b delayed by d (a[n] = b[n-d]), at k = +d.
    The cross spectrum X_a X_b* then has the phase +2 pi f d (b delayed)."""
    a = np.asarray(a); b = np.asarray(b)
    h = n_lags // 2
    lags = np.arange(-h, n_lags - h)
    n = min(a.size, b.size) - n_lags
    if n <= 0:
        raise ValueError("series shorter than the number of lags")
    bb = np.conj(b[h:h + n])
    R = np.array([np.mean(a[h + k:h + k + n] * bb) for k in lags])
    return lags, R


def xf_cross_spectrum(a: np.ndarray, b: np.ndarray, n_lags: int) -> np.ndarray:
    """XF correlator: Fourier transform of :func:`lag_correlation` (Wiener-
    Khinchin). n_lags channels in fftshift order (as ``Accumulator``); unit
    of the per-sample correlation, comparable to ``Accumulator.visibility()
    / window_power`` with the rectangular window."""
    lags, R = lag_correlation(a, b, n_lags)
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(R)))
