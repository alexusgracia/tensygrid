# build.ps1 — Compila la documentació Sphinx (equivalent a build.sh per a Windows)
#
# Ús (des de la carpeta manual/):
#   .\build.ps1
#
# Prerequisits:
#   - Python instal·lat i accessible des de PowerShell
#   - Entorn virtual creat: python -m venv venv
#   - Dependències instal·lades: pip install -r requirements.txt

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

# ---------------------------------------------------------------------------
# Activar entorn virtual si existeix
# ---------------------------------------------------------------------------
if (Test-Path "venv\Scripts\Activate.ps1") {
    . "venv\Scripts\Activate.ps1"
}

# ---------------------------------------------------------------------------
# Resolució del path al codi font
# ---------------------------------------------------------------------------
$CodePath = $env:TENSYGRID_CODE_PATH

if (-not $CodePath) {
    if (Test-Path "docs\local_paths.py") {
        $CodePath = python -c "import sys; sys.path.insert(0,'docs'); import local_paths; print(local_paths.TENSYGRID_CODE_PATH)"
    }
}

if (-not $CodePath) {
    Write-Warning "TENSYGRID_CODE_PATH no esta configurat."
    Write-Host "   La documentacio de l'API NO es generara."
    Write-Host "   Per activar-la:"
    Write-Host "     1) Copia docs\local_paths.py.example -> docs\local_paths.py i edita el path"
    Write-Host "     2) O defineix: `$env:TENSYGRID_CODE_PATH = 'C:\ruta\a\trunk\tensygrid'"
} else {
    Write-Host "OK  Codi font detectat a: $CodePath"
    $env:TENSYGRID_CODE_PATH = $CodePath
}

# ---------------------------------------------------------------------------
# Compilació
# ---------------------------------------------------------------------------
Write-Host ""
Write-Host "Compilant documentacio..."
sphinx-build -b html docs docs/_build/html

Write-Host ""
Write-Host "OK  Build completat: docs/_build/html/"
Write-Host "    Obre docs/_build/html/index.html al navegador per previsualitzar."
