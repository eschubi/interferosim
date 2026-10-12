# build_api_docs.ps1 -- generate the API documentation with Doxygen (Windows, PowerShell)
#
# Usage (in the project folder):
#     powershell -ExecutionPolicy Bypass -File build_api_docs.ps1
#     powershell -ExecutionPolicy Bypass -File build_api_docs.ps1 -Pdf -NoOpen
# or double-click build_api_docs.bat.
#
# Steps: find doxygen (PATH or the default installation folder); if it is
# missing, install it with winget (package DimitriVanHeesch.Doxygen); run
# "doxygen Doxyfile"; optionally build the PDF reference manual; open the
# HTML start page.
#
# Version 1.0.0 (2026-10-12)

param(
    [switch]$Pdf,        # also build docs\api\latex\refman.pdf (needs LaTeX)
    [switch]$NoOpen,     # do not open the HTML page in the browser
    [switch]$NoInstall   # do not install doxygen automatically
)

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

function Find-Doxygen {
    $cmd = Get-Command doxygen -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    foreach ($p in @("$env:ProgramFiles\doxygen\bin\doxygen.exe",
                     "${env:ProgramFiles(x86)}\doxygen\bin\doxygen.exe",
                     "$env:LOCALAPPDATA\Programs\doxygen\bin\doxygen.exe")) {
        if ($p -and (Test-Path $p)) { return $p }
    }
    return $null
}

# ---- 1. doxygen ----------------------------------------------------------
$dox = Find-Doxygen
if (-not $dox) {
    if ($NoInstall) { throw "doxygen not found; install it from https://www.doxygen.nl/download.html" }
    if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
        throw "doxygen not found and winget not available; install doxygen from https://www.doxygen.nl/download.html"
    }
    Write-Host "doxygen not found -- installing with winget (DimitriVanHeesch.Doxygen)"
    winget install --id DimitriVanHeesch.Doxygen -e --accept-package-agreements --accept-source-agreements
    if ($LASTEXITCODE -ne 0) { throw "winget install failed ($LASTEXITCODE)" }
    $dox = Find-Doxygen
    if (-not $dox) { throw "doxygen installed but not found; open a new PowerShell window and run the script again" }
}
Write-Host "doxygen: $dox  ($(& $dox --version))"

# ---- 2. generate -----------------------------------------------------------
New-Item -ItemType Directory -Force -Path "docs\api" | Out-Null
& $dox Doxyfile
if ($LASTEXITCODE -ne 0) { throw "doxygen failed ($LASTEXITCODE)" }
$warn = "docs\api\doxygen_warnings.txt"
if ((Test-Path $warn) -and ((Get-Item $warn).Length -gt 0)) {
    Write-Warning "doxygen reported warnings, see $warn"
}
Write-Host "HTML:  docs\api\html\index.html"

# ---- 3. optional PDF -----------------------------------------------------
if ($Pdf) {
    Push-Location "docs\api\latex"
    try {
        cmd /c make.bat
        if (Test-Path "refman.pdf") { Write-Host "PDF:   docs\api\latex\refman.pdf" }
        else { Write-Warning "refman.pdf was not created, see the LaTeX output above" }
    } finally { Pop-Location }
}

# ---- 4. open -------------------------------------------------------------
if (-not $NoOpen) { Start-Process "docs\api\html\index.html" }
