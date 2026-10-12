# src/interferosim/_log.py -- logger setup of interferosim
"""Central logger access. All modules obtain their logger via get_logger().

The root name is "interferosim" (aas6sim used "aas6").
"""

import logging

LOGGER_ROOT = "interferosim"


def get_logger(name: str) -> logging.Logger:
    """Return a child logger below the package root logger."""
    if name.startswith(LOGGER_ROOT):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_ROOT}.{name}")
