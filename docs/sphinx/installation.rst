.. docs/sphinx/installation.rst -- installation
.. Version 1.0.1 (2026-11-09): pip installation and uv development environment added
.. 2026-11-09: activation command per Windows shell (PowerShell Activate.ps1, cmd activate.bat)

Installation
============

Python 3.11 or newer is required.

Installation with pip
---------------------

The package installs from a fresh clone with plain ``pip``; ``uv`` is not
needed for this:

.. code-block:: text

   git clone https://github.com/eschubi/interferosim
   cd interferosim
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1               # Windows PowerShell
                                               # (cmd.exe: .venv\Scripts\activate.bat;
                                               #  Linux/macOS: source .venv/bin/activate)
   python -m pip install .                     # runtime dependencies only
   python -m pip install ".[notebooks]"        # optional: JupyterLab, ipywidgets

Development environment with uv
-------------------------------

The development environment is managed with `uv <https://docs.astral.sh/uv/>`_
and reproduced exactly from the lockfile ``uv.lock`` in the repository:

.. code-block:: text

   uv sync                  # .venv with runtime dependencies + group "dev"
   uv sync --group docs     # additionally the documentation tools
   uv run pytest
   uv run jupyter lab notebooks

Development dependencies are PEP 735 dependency groups
(``[dependency-groups]`` in ``pyproject.toml``): ``dev`` (JupyterLab,
ipywidgets, ipykernel and the group ``test``), ``test`` (pytest, nbval) and
``docs``. ``pip install .`` never installs them. After a change of
dependencies run ``uv lock`` and commit ``uv.lock`` with ``pyproject.toml``;
``.venv/`` is never committed.

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
