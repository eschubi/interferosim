@echo off
rem build_api_docs.bat -- generate the API documentation (Windows, double-click)
rem Calls build_api_docs.ps1 with all arguments, e.g.  build_api_docs.bat -Pdf
rem Version 1.0.0 (2026-10-12)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0build_api_docs.ps1" %*
if errorlevel 1 (
    echo.
    echo build_api_docs failed, see the messages above.
)
pause
