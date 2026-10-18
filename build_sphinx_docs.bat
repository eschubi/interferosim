@echo off
rem build_sphinx_docs.bat -- build the Sphinx documentation (Windows, double-click)
rem Output: docs\sphinx\_build\html\index.html
rem Version 1.0.0 (2026-10-18)
set "PY=%~dp0.venv\Scripts\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -c "import sphinx, furo, myst_parser, sphinx_copybutton" 2>nul || "%PY%" -m pip install -e "%~dp0.[docs]"
"%PY%" -m sphinx -b html "%~dp0docs\sphinx" "%~dp0docs\sphinx\_build\html"
if errorlevel 1 (
    echo.
    echo build_sphinx_docs failed, see the messages above.
) else (
    start "" "%~dp0docs\sphinx\_build\html\index.html"
)
pause
