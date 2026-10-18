.. docs/sphinx/notebooks.rst -- teaching notebooks

Notebooks
=========

The teaching notebooks in ``notebooks/`` use the package and add
interactive selectors (ipywidgets). Start them with

.. code-block:: text

   jupyter lab notebooks          # kernel "Python (interferosim)"

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Notebook
     - Content
   * - ``01_sky_models.ipynb``
     - sky models, file format, maps, own skies, visibilities, catalogue of
       known regions with selector
   * - ``02_instrument.ipynb``
     - antenna types and primary beams, dish transfer function, B210
       receiver and sensitivity, arrays and baselines, choice of the
       baseline length (Cyg A), selector

The migration of Appendix A (Wiener filter vs. CLEAN) and Appendix B (the
Sun as a first target) follows in the next steps. The ``binder/`` folder
holds the environment for Binder / JupyterHub.
