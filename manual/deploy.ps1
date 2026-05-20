# deploy.ps1 — Compila la documentació i la desplegua a GitHub Pages (equivalent a deploy.sh per a Windows)
#
# Ús (des de la carpeta manual/):
#   .\deploy.ps1              # compila + puja a gh-pages
#   .\deploy.ps1 -BuildOnly   # només compila, sense pujar
#
# Prerequisits:
#   - Python instal·lat i accessible des de PowerShell
#   - Entorn virtual creat: python -m venv venv
#   - Dependències instal·lades: pip install -r requirements.txt  (inclou ghp-import)
#   - Accés de push al repositori remot origin

param(
    [switch]$BuildOnly
)

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
    Write-Host "   La seccio de l'API NO es generara."
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

# Fitxer .nojekyll per evitar que GitHub Pages processi els estils amb Jekyll
New-Item -ItemType File -Path "docs\_build\html\.nojekyll" -Force | Out-Null

Write-Host "OK  Build completat: docs/_build/html/"

# ---------------------------------------------------------------------------
# Deploy a gh-pages
# ---------------------------------------------------------------------------
if ($BuildOnly) {
    Write-Host ""
    Write-Host "Mode -BuildOnly: s'ha omes el desplegament."
    exit 0
}

if (-not (Get-Command "ghp-import" -ErrorAction SilentlyContinue)) {
    Write-Error "ghp-import no trobat. Instal·la'l amb: pip install ghp-import"
    exit 1
}

Write-Host ""
Write-Host "Desplegant a branca gh-pages..."
ghp-import -n -p -f docs\_build\html

Write-Host ""
Write-Host "OK  Documentacio publicada a GitHub Pages!"
Write-Host "   La pagina pot trigar uns instants a actualitzar-se."
