@echo off
rem setup_venv.bat -- create the project virtual environment (Windows, double-click)
rem Calls setup_venv.ps1 with all arguments, e.g.  setup_venv.bat -Recreate
rem Version 1.0.0 (2026-10-12)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup_venv.ps1" %*
if errorlevel 1 (
    echo.
    echo setup_venv failed, see the messages above.
)
pause
