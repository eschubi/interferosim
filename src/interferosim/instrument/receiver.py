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

Version 0.3.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, Optional, Union

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
        """Read ``configs/receivers/<name>.toml``."""
        with open(path, "rb") as fh:
            return cls.from_dict(tomllib.load(fh))
