#!/bin/bash
# ---------------------------------------------------------------------------
# deploy.sh — Compila la documentació i la desplegua a GitHub Pages
#
# Flux de treball:
#   1. Comprova que TENSYGRID_CODE_PATH és accessible (per generar l'API)
#   2. Compila amb sphinx-build → manual/docs/_build/html/
#   3. Puja automàticament el HTML a la branca `gh-pages` del remot origin
#
# Prerequisits:
#   - pip install ghp-import   (o: pip install -r requirements.txt)
#   - Tenir acces de push al repositori remot origin
#
# Ús:
#   cd manual
#   ./deploy.sh                 # compila + desplegua
#   ./deploy.sh --build-only    # només compila, sense pujar
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BUILD_ONLY=false
if [[ "$1" == "--build-only" ]]; then
    BUILD_ONLY=true
fi

# ---------------------------------------------------------------------------
# 1. Resolució del path al codi font
# ---------------------------------------------------------------------------
if [ -z "$TENSYGRID_CODE_PATH" ]; then
    if [ -f "docs/local_paths.py" ]; then
        TENSYGRID_CODE_PATH=$(python3 -c \
            "import sys; sys.path.insert(0,'docs'); import local_paths; print(local_paths.TENSYGRID_CODE_PATH)")
    fi
fi

if [ -z "$TENSYGRID_CODE_PATH" ]; then
    echo "⚠️  AVÍS: TENSYGRID_CODE_PATH no està configurat."
    echo "   La secció de l'API NO es generarà."
    echo "   Per activar-la:"
    echo "     1) Copia docs/local_paths.py.example → docs/local_paths.py i edita el path"
    echo "     2) O exporta: export TENSYGRID_CODE_PATH=/ruta/a/trunk/tensygrid"
    echo ""
else
    echo "✅  Codi font detectat a: $TENSYGRID_CODE_PATH"
    export TENSYGRID_CODE_PATH
fi

# ---------------------------------------------------------------------------
# 2. Compilació
# ---------------------------------------------------------------------------
echo ""
echo "🔨 Compilant documentació..."
sphinx-build -b html docs docs/_build/html
touch docs/_build/html/.nojekyll    # evita que Jekyll processi els fitxers

echo "✅  Build completat: docs/_build/html/"

# ---------------------------------------------------------------------------
# 3. Deploy a gh-pages
# ---------------------------------------------------------------------------
if [ "$BUILD_ONLY" = true ]; then
    echo ""
    echo "ℹ️  Mode --build-only: s'ha omès el desplegament."
    exit 0
fi

if ! command -v ghp-import &> /dev/null; then
    echo ""
    echo "❌  ghp-import no trobat. Instal·la'l amb:"
    echo "      pip install ghp-import"
    echo "    o:"
    echo "      pip install -r requirements.txt"
    exit 1
fi

echo ""
echo "🚀 Desplegant a branca gh-pages..."
# -n : no-jekyll (.nojekyll ja creat manualment)
# -p : push automàtic a origin/gh-pages
# -f : força (substitueix commits anteriors de gh-pages)
ghp-import -n -p -f docs/_build/html

echo ""
echo "✅  Documentació publicada a GitHub Pages!"
echo "   La pàgina pot trigar uns instants a actualitzar-se."
