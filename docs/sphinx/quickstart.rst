.. docs/sphinx/quickstart.rst -- quick start

Quick start
===========

Activate the environment, draw a sky model and open the notebooks:

.. code-block:: text

   .venv\Scripts\activate          # Windows   (Linux/macOS: source .venv/bin/activate)
   interferosim-sky configs/skies/aas6_demo.json --overview -o demo.png
   jupyter lab notebooks          # kernel "Python (interferosim)"

Command line
------------

``interferosim-sky`` shows or saves the map of a sky model file
(:mod:`interferosim.sky.cli`):

.. code-block:: text

   interferosim-sky configs/skies/aas6_demo.json                 # window
   interferosim-sky configs/skies/aas6_demo.json -o demo.png     # file
   interferosim-sky configs/skies/aas6_demo.json --overview -o demo.png
   interferosim-sky --demo -o demo.png                           # built-in demo

From Python
-----------

.. code-block:: python

   from interferosim.sky import SkyModel
   from interferosim.instrument import Antenna, Receiver, Array

   sky = SkyModel.from_json("configs/skies/cyg_a.json")
   img, x, y, pix = sky.render_map(npix=512)        # Jy per pixel

   dish = Antenna.from_toml("configs/antennas/dish_3m.toml")
   b210 = Receiver.from_toml("configs/receivers/b210.toml")
   arr  = Array.from_toml("configs/sites/aas6_sim_stations.toml")

   b210.sefd_jy(dish)                               # SEFD in Jy
   arr.baselines_enu()                              # baselines east/north/up in m
   b210.baseline_sigma_jy(dish, dish, 10.0)         # noise per baseline, 10 s
