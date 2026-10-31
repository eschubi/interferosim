# src/interferosim/signal/fx.py -- FX correlator (interferosim.signal.fx)
"""FX correlator for the IQ streams of all stations.

"F" -- every stream is cut into windows of ``nfft`` samples and Fourier
transformed; "X" -- the spectra of each pair of stations are multiplied,
X_i X_j*, and averaged over the integration time. The result is a cross
spectrum per baseline and integration: a :class:`~interferosim.visdata.VisData`
set with ``nfft`` frequency channels.

Processing chain (as AAS6 ``tools/iq_correlate.py`` and
``modules/interfero/correlate.py``):

1. **Model correction** (delay tracking and fringe stopping with the
   geometric model of the phase centre), :func:`correct_geometry`.
2. **Residual fringe solution** per station against the reference --
   clock offset (residual delay) and LO offset (residual fringe rate) are
   unknown to the model; :meth:`FXCorrelator.search` finds them with
   :func:`delay_rate_search` on a bright segment,
   :func:`apply_fringe_solution` removes them.
3. **Cross spectra** per integration with :class:`Accumulator`; normalised
   by the window power, so that the visibility is in Jy per channel, with
   the weight 1/sigma^2 from the auto spectra:
   ``sigma^2 = <|X_i|^2> <|X_j|^2> / (2 n_win)`` (real and imaginary part each).

Version 0.6.0 (2026-10-31): first implementation.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from ..observation.geometry import SIDEREAL_DEG_PER_S
from ..visdata import VisData
from .correlate import (Accumulator, FringeSolution, apply_fringe_solution, correct_geometry,
                        delay_rate_search, refine_fringe_solution)

__all__ = ["FXCorrelator"]


class FXCorrelator:
    """FX correlator with ``nfft`` channels and integration time ``t_int_s``."""

    def __init__(self, nfft: int = 64, t_int_s: float = 1.0, window: str = "hann"):
        if nfft < 2:
            raise ValueError("nfft must be >= 2")
        self.nfft = int(nfft)
        self.t_int_s = float(t_int_s)
        self.window = window

    # ------------------------------------------------------------------
    def correct(self, blocks: Dict[str, "object"], sim, lha0_deg: float = 0.0,
                solutions: Optional[Dict[str, FringeSolution]] = None,
                model: bool = True) -> Dict[str, np.ndarray]:
        """Streams of all stations aligned to the reference: model delay and
        fringe phase removed (``model``), then the residual fringe solution
        of each station (``solutions`` by station name)."""
        solutions = solutions or {}
        out = {}
        for k, st in enumerate(sim.array.stations):
            x = np.asarray(blocks[st.name].data)
            if model and k > 0:
                x = correct_geometry(x, lambda t, k=k: sim.phase_centre_delay_s(k, t, lha0_deg),
                                     sim.fs, sim.fc)
            if st.name in solutions:
                x = apply_fringe_solution(x, solutions[st.name], sim.fs)
            out[st.name] = x
        return out

    def search(self, blocks: Dict[str, "object"], sim, station: str, *, segment_s: float = 2.0,
               start_s: float = 0.0, lha0_deg: float = 0.0, model: bool = True,
               max_delay_samples: float = 64.0, max_rate_hz: float = 5.0,
               n_rate: int = 201, refine: bool = True) -> FringeSolution:
        """Residual delay and fringe rate of ``station`` against the reference
        on the segment [start_s, start_s + segment_s] after the model correction;
        with ``refine`` the solution is refined by phase fits
        (:func:`refine_fringe_solution`)."""
        ref = sim.array.stations[0].name
        streams = self.correct(blocks, sim, lha0_deg, model=model)
        n0 = int(start_s * sim.fs)
        n = int(segment_s * sim.fs)
        a, b = streams[ref][n0:n0 + n], streams[station][n0:n0 + n]
        sol = delay_rate_search(a, b, sim.fs, max_delay_samples=max_delay_samples,
                                max_rate_hz=max_rate_hz, n_rate=n_rate)
        if refine:
            sol = refine_fringe_solution(a, b, sol, sim.fs, nfft=self.nfft)
        return sol

    # ------------------------------------------------------------------
    def correlate(self, blocks: Dict[str, "object"], sim, lha0_deg: float = 0.0,
                  solutions: Optional[Dict[str, FringeSolution]] = None, model: bool = True,
                  start_mjd: float = 0.0) -> VisData:
        """Cross spectra of all baselines, one row per integration and
        baseline, ``nfft`` channels."""
        streams = self.correct(blocks, sim, lha0_deg, solutions, model)
        st = sim.array.stations
        n_total = min(len(v) for v in streams.values())
        n_int = max(int(round(self.t_int_s * sim.fs)) // self.nfft, 1) * self.nfft
        n_rows_t = n_total // n_int
        if n_rows_t == 0:
            raise ValueError("recording shorter than one integration")
        acc0 = Accumulator.create(self.nfft, self.window)
        freq = acc0.freq_axis(sim.fs, sim.fc)
        wpow = acc0.window_power
        rows = {k: [] for k in ("t", "ha", "u", "v", "w", "b", "vis", "wt")}
        t_mid = (np.arange(n_rows_t) + 0.5) * n_int / sim.fs
        for b_idx, (i, j) in enumerate(sim.array.baseline_pairs()):
            a_s, b_s = streams[st[i].name], streams[st[j].name]
            q = sim.uvw(t_mid, lha0_deg)[b_idx] * (freq[0] / sim.fc)   # at freq_hz[0]
            for r in range(n_rows_t):
                acc = Accumulator.create(self.nfft, self.window)
                sl = slice(r * n_int, (r + 1) * n_int)
                acc.add(a_s[sl], b_s[sl])
                vis = acc.visibility() / wpow
                pa = acc.auto_a / acc.n_win / wpow
                pb = acc.auto_b / acc.n_win / wpow
                var = pa * pb / (2.0 * acc.n_win)
                rows["vis"].append(vis); rows["wt"].append(1.0 / np.maximum(var, 1e-300))
            rows["t"].append(t_mid); rows["ha"].append(lha0_deg + SIDEREAL_DEG_PER_S * t_mid)
            rows["u"].append(q[:, 0]); rows["v"].append(q[:, 1]); rows["w"].append(q[:, 2])
            rows["b"].append(np.full(n_rows_t, b_idx))
        cat = lambda k: np.concatenate(rows[k])
        return VisData(time_mjd=start_mjd + cat("t") / 86400.0, ha_deg=cat("ha"),
                       u=cat("u"), v=cat("v"), w=cat("w"), baseline=cat("b").astype(int),
                       vis=np.array(rows["vis"]), weight=np.array(rows["wt"]), freq_hz=freq,
                       baseline_names=sim.array.baseline_names(),
                       header={"OBJECT": sim.sky.center.label or sim.sky.name,
                               "RA": sim.sky.center.ra_deg, "DEC": sim.sky.center.dec_deg,
                               "TELESCOP": sim.array.name, "INSTRUME": sim.receiver.name,
                               "NFFT": self.nfft, "INTTIME": n_int / sim.fs, "WINDOW": str(self.window),
                               "MODELCOR": bool(model), "FRINGSOL": sorted((solutions or {}).keys()),
                               "SIMULATD": True, "ORIGIN": "interferosim.signal.FXCorrelator"})
