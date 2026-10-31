# src/interferosim/instrument/receiver.py -- receiver model (interferosim.instrument.receiver)
"""Receiver model with the device profile of the Ettus USRP B210.

The receiver provides the system temperature, the bandwidth and the
sampling parameters; together with an antenna it gives the system
equivalent flux density (SEFD) and the radiometer noise.

Origin
------
``aas6sim`` 1.3.1 ``model.Instrument`` (``t_sys_k``, ``bandwidth_hz``,
``sefd_jy``) and ``simulate.py`` (noise per baseline
sigma = SEFD / sqrt(2 B tau)); device settings of notebook B
(``configs/receivers/b210.toml``).

Formulas
--------
* SEFD = 2 k T_sys / A_eff                                  [Jy]
* single dish (one polarisation): sigma_S = SEFD / sqrt(B tau)
* baseline i-j, real and imaginary part each:
  sigma = sqrt(SEFD_i SEFD_j) / sqrt(2 B tau)   (aas6sim convention)

Device profile (since 0.3.2)
----------------------------
Besides the operating settings a profile may carry the data-sheet values
of the device (frequency range, ADC resolution, maximum sample rate and
bandwidth, noise figure, gain range, frequency accuracy of the reference,
PPS accuracy). They give

* the receiver noise temperature T_rx = 290 K (10^(NF/10) - 1),
* the data rate of the recorded IQ stream,
* the LO offset y f and the time t_1 = 1 / (2 pi y f) after which the
  phase has drifted by one radian (relative frequency error y of the
  reference, chosen by ``clock_source``),
* :meth:`Receiver.validate` -- warnings when a setting exceeds the device.

A profile may name another profile with ``base = "<name>"`` (same
directory); its keys are read first and overridden by the file
(e.g. ``b210_gpsdo.toml`` builds on ``b210.toml``).

Version 0.3.2 (2026-10-31): device profile (data-sheet fields), data rate,
  LO offset / coherence time, receiver noise temperature, validate(),
  profiles with ``base``.
Version 0.3.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ..constants import HI_21CM, JANSKY, K_BOLTZMANN
from .antennas import Antenna

__all__ = ["Receiver"]


@dataclass
class Receiver:
    """Receiver (SDR) with system temperature and sampling parameters."""

    name: str = "receiver"
    t_sys_k: float = 120.0
    bandwidth_hz: float = 2.0e6
    sample_rate_sps: float = 2.0e6
    rx_channels: int = 1
    bits: int = 16
    agc: bool = False
    dc_offset_correction: bool = False
    iq_balance_correction: bool = False
    clock_source: str = "internal"
    # ---- operating settings added in 0.3.2 (None = not set) ----------
    center_freq_hz: Optional[float] = None      # tuning frequency
    rx_gain_db: Optional[float] = None          # manual RX gain
    sample_format: str = "sc16"                 # host format: sc16, sc8, fc32
    time_source: str = "none"                   # none, gpsdo, external (PPS)
    external_lna: bool = False                  # LNA in front of the receiver (t_sys_k includes it)
    # ---- device data (data sheet), added in 0.3.2 ---------------------
    freq_min_hz: Optional[float] = None
    freq_max_hz: Optional[float] = None
    adc_bits: Optional[int] = None
    max_sample_rate_sps: Optional[float] = None         # one channel
    max_sample_rate_2ch_sps: Optional[float] = None     # two channels (MIMO)
    max_bandwidth_hz: Optional[float] = None            # one channel
    max_bandwidth_2ch_hz: Optional[float] = None        # two channels
    noise_figure_db: Optional[float] = None
    rx_gain_range_db: Optional[List[float]] = None
    iip3_dbm: Optional[float] = None
    adc_sfdr_dbc: Optional[float] = None
    freq_accuracy_ppm: Optional[float] = None           # internal reference
    freq_accuracy_gpsdo_unlocked_ppb: Optional[float] = None
    freq_accuracy_gpsdo_locked_ppb: Optional[float] = None
    pps_accuracy_ns: Optional[float] = None             # locked GPSDO
    host_interface: Optional[str] = None
    fpga: Optional[str] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    def sefd_jy(self, antenna: Antenna, f_hz: float = HI_21CM) -> float:
        """System equivalent flux density of ``antenna`` with this receiver."""
        return 2.0 * K_BOLTZMANN * self.t_sys_k / antenna.effective_area(f_hz) / JANSKY

    def single_dish_sigma_jy(self, antenna: Antenna, t_int_s: float,
                             f_hz: float = HI_21CM,
                             bandwidth_hz: Optional[float] = None) -> float:
        """Radiometer noise of one antenna, SEFD / sqrt(B tau)."""
        b = bandwidth_hz or self.bandwidth_hz
        return self.sefd_jy(antenna, f_hz) / math.sqrt(b * t_int_s)

    def baseline_sigma_jy(self, antenna_1: Antenna, antenna_2: Antenna, t_int_s: float,
                          f_hz: float = HI_21CM,
                          bandwidth_hz: Optional[float] = None) -> float:
        """Noise of the real (and of the imaginary) part of one visibility,
        sqrt(SEFD_1 SEFD_2) / sqrt(2 B tau)."""
        b = bandwidth_hz or self.bandwidth_hz
        s = math.sqrt(self.sefd_jy(antenna_1, f_hz) * self.sefd_jy(antenna_2, f_hz))
        return s / math.sqrt(2.0 * b * t_int_s)

    @property
    def channel_width_hz(self) -> float:
        """Sample rate of one complex channel = width of the analysed band."""
        return self.sample_rate_sps

    # ------------------------------------------------------------------
    # device profile (0.3.2)
    _BYTES_PER_SAMPLE = {"sc8": 2, "sc16": 4, "fc32": 8}

    @property
    def bytes_per_sample(self) -> int:
        """Bytes of one complex sample in the host format."""
        try:
            return self._BYTES_PER_SAMPLE[self.sample_format]
        except KeyError:
            raise ValueError(f"unknown sample_format {self.sample_format!r}") from None

    def data_rate_bytes_per_s(self, channels: Optional[int] = None) -> float:
        """Data rate of the recorded IQ stream (default: all RX channels)."""
        n = self.rx_channels if channels is None else channels
        return self.sample_rate_sps * self.bytes_per_sample * n

    def receiver_noise_temperature_k(self) -> Optional[float]:
        """Noise temperature of the receiver from its noise figure,
        T_rx = 290 K (10^(NF/10) - 1); None without noise figure."""
        if self.noise_figure_db is None:
            return None
        return 290.0 * (10.0 ** (self.noise_figure_db / 10.0) - 1.0)

    def relative_frequency_error(self, clock_source: Optional[str] = None) -> Optional[float]:
        """Relative frequency error y of the reference for ``clock_source``
        ("internal", "gpsdo", "gpsdo_unlocked"); None if not known."""
        src = clock_source or self.clock_source
        if src == "internal" and self.freq_accuracy_ppm is not None:
            return self.freq_accuracy_ppm * 1e-6
        if src == "gpsdo" and self.freq_accuracy_gpsdo_locked_ppb is not None:
            return self.freq_accuracy_gpsdo_locked_ppb * 1e-9
        if src == "gpsdo_unlocked" and self.freq_accuracy_gpsdo_unlocked_ppb is not None:
            return self.freq_accuracy_gpsdo_unlocked_ppb * 1e-9
        return None

    def lo_offset_hz(self, f_hz: Optional[float] = None,
                     clock_source: Optional[str] = None) -> Optional[float]:
        """Largest LO offset y f at ``f_hz`` (default: centre frequency)."""
        y = self.relative_frequency_error(clock_source)
        f = f_hz or self.center_freq_hz or HI_21CM
        return None if y is None else y * f

    def phase_coherence_time_s(self, f_hz: Optional[float] = None,
                               clock_source: Optional[str] = None) -> Optional[float]:
        """Time after which an LO offset y f has turned the phase by one
        radian, t_1 = 1 / (2 pi y f). A lower limit: the data-sheet value is
        an accuracy, the short-term stability is usually better."""
        df = self.lo_offset_hz(f_hz, clock_source)
        return None if not df else 1.0 / (2.0 * math.pi * df)

    def validate(self) -> List[str]:
        """Warnings where the settings exceed the device data (empty list:
        everything consistent). Data not given in the profile is not checked."""
        w: List[str] = []
        two = self.rx_channels >= 2
        max_sr = self.max_sample_rate_2ch_sps if two else self.max_sample_rate_sps
        if max_sr and self.sample_rate_sps > max_sr:
            w.append(f"sample rate {self.sample_rate_sps:g} S/s exceeds the maximum "
                     f"{max_sr:g} S/s for {self.rx_channels} channel(s)")
        max_bw = self.max_bandwidth_2ch_hz if two else self.max_bandwidth_hz
        if max_bw and self.bandwidth_hz > max_bw:
            w.append(f"bandwidth {self.bandwidth_hz:g} Hz exceeds the maximum {max_bw:g} Hz")
        if self.bandwidth_hz > self.sample_rate_sps:
            w.append("bandwidth larger than the complex sample rate (aliasing)")
        f = self.center_freq_hz
        if f is not None and self.freq_min_hz and self.freq_max_hz \
                and not (self.freq_min_hz <= f <= self.freq_max_hz):
            w.append(f"centre frequency {f:g} Hz outside {self.freq_min_hz:g}-{self.freq_max_hz:g} Hz")
        g, rng = self.rx_gain_db, self.rx_gain_range_db
        if g is not None and rng and not (rng[0] <= g <= rng[1]):
            w.append(f"RX gain {g:g} dB outside the range {rng[0]:g}-{rng[1]:g} dB")
        t_rx = self.receiver_noise_temperature_k()
        if t_rx is not None and self.t_sys_k < t_rx and not self.external_lna:
            w.append(f"t_sys_k = {self.t_sys_k:g} K is below the receiver noise temperature "
                     f"{t_rx:.0f} K (NF {self.noise_figure_db:g} dB): only possible with an "
                     f"external LNA (set external_lna = true)")
        if self.clock_source == "internal" and self.time_source != "none":
            w.append("time source set but internal clock: the frequency reference is free-running")
        return w

    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        d = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "extra"}
        d.update(self.extra)
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Receiver":
        allowed = {f.name for f in fields(cls)} - {"extra"}
        kwargs = {k: v for k, v in d.items() if k in allowed}
        extra = {k: v for k, v in d.items() if k not in allowed}
        return cls(**kwargs, extra=extra)

    @classmethod
    def from_toml(cls, path: Union[str, Path]) -> "Receiver":
        """Read ``configs/receivers/<name>.toml``.

        A key ``base = "<name>"`` reads ``<name>.toml`` of the same directory
        first; the keys of the file override it (since 0.3.2).
        """
        return cls.from_dict(cls._load_toml(Path(path)))

    @classmethod
    def _load_toml(cls, path: Path, _seen: Optional[set] = None) -> Dict[str, Any]:
        seen = _seen or set()
        if path.resolve() in seen:
            raise ValueError(f"circular 'base' in {path}")
        seen.add(path.resolve())
        with open(path, "rb") as fh:
            d = tomllib.load(fh)
        base = d.pop("base", None)
        if base:
            merged = cls._load_toml(path.parent / f"{base}.toml", seen)
            merged.update(d)
            d = merged
        return d
