.. docs/sphinx/installation.rst -- installation

Installation
============

Python 3.11 or newer is required.

With the scripts
----------------

Create the project environment ``.venv``:

.. code-block:: text

   setup_venv.bat                 # Windows (double-click or command line)
   ./setup_venv.sh                # Linux / macOS

Development installation
------------------------

.. code-block:: text

   python -m venv .venv
   .venv\Scripts\python -m pip install -e ".[notebooks,test]"
   .venv\Scripts\python -m pytest

Optional dependency groups (``pyproject.toml``):

``notebooks``
   JupyterLab, ipywidgets, ipykernel
``test``
   pytest, nbval
``docs``
   Sphinx, the Furo theme, MyST parser, sphinx-copybutton

Building this documentation
---------------------------

.. code-block:: text

   .venv\Scripts\python -m pip install -e ".[docs]"
   build_sphinx_docs.bat                    # -> docs/sphinx/_build/html/index.html

or directly:

.. code-block:: text

   python -m sphinx -b html docs/sphinx docs/sphinx/_build/html

The Doxygen API documentation (``build_api_docs.bat``) and the LaTeX
documentation (``cd docs && latexmk``) remain available.
