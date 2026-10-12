# src/interferosim/sky/__init__.py -- sub-package interferosim.sky
"""Sky models: point, Gaussian, disc and limb-brightened disc components, HI and OH maser sources, rendering to images.

Version 0.2.0 (2026-10-12): components, SkyModel and the JSON format of
AAS6 ``json/interfero/sim_sky.json`` implemented (was: empty skeleton).

Example::

    from interferosim.sky import SkyModel
    sky = SkyModel.from_json("configs/skies/aas6_demo.json")
    img, x, y, pix = sky.render_map(npix=512)
"""

from .components import (COMPONENT_KINDS, DiscSource, GaussianSource,
                         PointSource, SkyComponent, component_from_dict)
from .models import AM, AS, FieldCenter, SkyModel

__all__ = [
    "SkyComponent", "PointSource", "GaussianSource", "DiscSource",
    "COMPONENT_KINDS", "component_from_dict",
    "FieldCenter", "SkyModel", "AS", "AM",
]
