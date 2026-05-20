#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# Comprovació del path al codi font de tensygrid (projecte A)
# ---------------------------------------------------------------------------
if [ -z "$TENSYGRID_CODE_PATH" ]; then
    # Intenta llegir-lo del fitxer local_paths.py si existeix
    if [ -f "docs/local_paths.py" ]; then
        TENSYGRID_CODE_PATH=$(python3 -c "import sys; sys.path.insert(0,'docs'); import local_paths; print(local_paths.TENSYGRID_CODE_PATH)")
    fi
fi

if [ -z "$TENSYGRID_CODE_PATH" ]; then
    echo "⚠️  AVÍS: TENSYGRID_CODE_PATH no està configurat."
    echo "   La documentació de l'API NO es generarà."
    echo "   Per activar-la:"
    echo "     1) Copia docs/local_paths.py.example → docs/local_paths.py i edita el path"
    echo "     2) O exporta: export TENSYGRID_CODE_PATH=/ruta/a/trunk/tensygrid"
else
    echo "✅  Codi font detectat a: $TENSYGRID_CODE_PATH"
    export TENSYGRID_CODE_PATH
fi

sphinx-build -b html docs docs/_build/html
