# setup_venv.ps1 -- create the project virtual environment .venv (Windows, PowerShell)
#
# Usage (in the project folder):
#     powershell -ExecutionPolicy Bypass -File setup_venv.ps1
#     powershell -ExecutionPolicy Bypass -File setup_venv.ps1 -Python "C:\Path\to\python.exe"
#     powershell -ExecutionPolicy Bypass -File setup_venv.ps1 -Recreate -DropboxIgnore
# or simply double-click setup_venv.bat.
#
# Steps: find Python >= 3.11, create .venv, install interferosim in
# editable mode with the extras "notebooks" and "test", register the
# Jupyter kernel "Python (interferosim)", run the tests.
#
# Version 1.0.0 (2026-10-12)

param(
    [string]$Python = "",        # explicit interpreter; default: search
    [switch]$Recreate,           # delete an existing .venv first
    [switch]$NoTests,            # skip pytest at the end
    [switch]$DropboxIgnore       # mark .venv as ignored by Dropbox sync
)

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

function Invoke-Checked {
    # run a native command and stop on a non-zero exit code
    param([string]$Exe, [string[]]$Arguments)
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "command failed ($LASTEXITCODE): $Exe $($Arguments -join ' ')"
    }
}

function Test-Python {
    # returns $true if the interpreter exists and is Python >= 3.11
    param([string]$Exe, [string[]]$Pre)
    try {
        $out = & $Exe @Pre -c "import sys; print(sys.version_info >= (3, 11))" 2>$null
        return ($LASTEXITCODE -eq 0 -and "$out".Trim() -eq "True")
    } catch {
        return $false
    }
}

# ---- 1. find the interpreter ---------------------------------------
$exe = $null; $pre = @()
if ($Python) {
    if (-not (Test-Python $Python @())) { throw "$Python is not Python >= 3.11" }
    $exe = $Python
} else {
    $candidates = @(
        @("py", "-3.13"), @("py", "-3.12"), @("py", "-3.11"), @("py", "-3"),
        @("python"), @("python3")
    )
    foreach ($c in $candidates) {
        $cExe = $c[0]
        $cPre = @()
        if ($c.Count -gt 1) { $cPre = $c[1..($c.Count - 1)] }
        if (Get-Command $cExe -ErrorAction SilentlyContinue) {
            if (Test-Python $cExe $cPre) { $exe = $cExe; $pre = $cPre; break }
        }
    }
    if (-not $exe) { throw "no Python >= 3.11 found; install it or pass -Python <path>" }
}
$ver = & $exe @pre -c "import sys; print(sys.version.split()[0])"
Write-Host "Python: $exe $($pre -join ' ')  (version $ver)"

# ---- 2. create .venv -------------------------------------------------
if ($Recreate -and (Test-Path ".venv")) {
    Write-Host "removing existing .venv"
    Remove-Item -Recurse -Force ".venv"
}
if (-not (Test-Path ".venv\Scripts\python.exe")) {
    Write-Host "creating .venv"
    Invoke-Checked $exe ($pre + @("-m", "venv", ".venv"))
} else {
    Write-Host ".venv exists, updating it"
}
$vpy = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"

if ($DropboxIgnore) {
    # Dropbox does not sync folders that carry the stream com.dropbox.ignored
    try {
        Set-Content -Path ".venv" -Stream "com.dropbox.ignored" -Value 1
        Write-Host ".venv marked as ignored by Dropbox"
    } catch {
        Write-Warning "could not mark .venv for Dropbox: $_"
    }
}

# ---- 3. install ------------------------------------------------------
Invoke-Checked $vpy @("-m", "pip", "install", "--upgrade", "pip")
Invoke-Checked $vpy @("-m", "pip", "install", "-e", ".[notebooks,test]")

# ---- 4. Jupyter kernel -----------------------------------------------
Invoke-Checked $vpy @("-m", "ipykernel", "install", "--user",
                      "--name", "interferosim", "--display-name", "Python (interferosim)")

# ---- 5. tests -------------------------------------------------------
if (-not $NoTests) {
    Invoke-Checked $vpy @("-m", "pytest", "-q")
}

Write-Host ""
Write-Host "done. Activate with:   .venv\Scripts\Activate.ps1"
Write-Host "Notebooks:             .venv\Scripts\jupyter-lab.exe notebooks"
Write-Host "Sky map:               .venv\Scripts\interferosim-sky.exe --demo --overview"
