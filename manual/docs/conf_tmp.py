import sys
import os
import warnings
sys.path.insert(0, os.path.abspath('.'))

# ---------------------------------------------------------------------------
# Path al codi font del projecte A (tensygrid)
#
# Cada usuari ha de configurar el path de la seva màquina amb UNA d'aquestes
# dues opcions (per ordre de prioritat):
#
#   1) Variable d'entorn:   export TENSYGRID_CODE_PATH=/ruta/a/trunk/tensygrid
#      (Windows):           set TENSYGRID_CODE_PATH=C:\ruta\a\trunk\tensygrid
#
#   2) Fitxer local (gitignored): copia local_paths.py.example → local_paths.py
#      i edita el path.
# ---------------------------------------------------------------------------
_code_path = os.environ.get('TENSYGRID_CODE_PATH', '')

if not _code_path:
    try:
        import local_paths
        _code_path = local_paths.TENSYGRID_CODE_PATH
    except ImportError:
        pass

if not _code_path:
    warnings.warn(
        "\n\n[conf.py] TENSYGRID_CODE_PATH no està configurat.\n"
        "La documentació de l'API no es generarà.\n"
        "Copia local_paths.py.example → local_paths.py i edita el path,\n"
        "o exporta la variable d'entorn TENSYGRID_CODE_PATH.\n"
    )

# ---------------------------------------------------------------------------
# Informació del projecte
# ---------------------------------------------------------------------------
project = 'Manuals'
author = 'Carlos Herrera Vázquez, Alexandre Gràcia-Calvo & Eduardo Prieto-Araujo'
copyright = '2026, Carlos Herrera Vázquez, Alexandre Gràcia-Calvo & Eduardo Prieto-Araujo'
release = '0.1'
version = '0.1'

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------
extensions = [
    'myst_parser',
    'sphinxcontrib.mermaid',
]

if _code_path:
    extensions += [
        'sphinx.ext.napoleon',      # Ha d'anar ABANS d'autoapi per processar docstrings NumPy/Google
        'sphinx.ext.viewcode',      # Enllaç al codi font
        'autoapi.extension',        # Documentació automàtica del codi Python
    ]

# ---------------------------------------------------------------------------
# Configuració de sphinx-autoapi
# ---------------------------------------------------------------------------
autoapi_type = 'python'
autoapi_dirs = [_code_path] if _code_path else []
autoapi_root = 'api'
autoapi_add_toctree_entry = True
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]
autoapi_ignore = [
    '*/backups/*',
    '*/__pycache__/*',
    '*/files/*',
    '*/precomputed_builds/*',
]
autoapi_python_class_content = 'both'   # docstring classe + __init__
autoapi_keep_files = True              # regenera sempre en cada build
autoapi_template_dir = '_templates'    # plantilles personalitzades amb filtre napoleon

# Evitar que Sphinx intenti parsejar les plantilles Jinja com a documents RST
exclude_patterns = ['_templates', '_build']

# Registrar el filtre Jinja que transforma docstrings NumPy/Google → RST via napoleon
def autoapi_prepare_jinja_env(jinja_env):
    from sphinx.ext.napoleon import Config
    from sphinx.ext.napoleon.docstring import NumpyDocstring, GoogleDocstring

    napoleon_config = Config(
        napoleon_numpy_docstring=True,
        napoleon_google_docstring=True,
        napoleon_use_param=True,
        napoleon_use_rtype=True,
        napoleon_use_ivar=True,
    )

    def napoleon_filter(docstring):
        if not docstring:
            return ''
        try:
            return str(NumpyDocstring(docstring, napoleon_config))
        except Exception:
            return docstring

    jinja_env.filters['napoleon'] = napoleon_filter

# Napoleon: format dels docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# MyST: activa extensions útils (taules, notes, enllaços automàtics...)
myst_enable_extensions = [
    'colon_fence',
    'deflist',
    'tasklist',
    'smartquotes',
    'strikethrough',
]

# ---------------------------------------------------------------------------
# Tema i opcions visuals
# ---------------------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    # Navegació
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    # Amplada de la pàgina (s'adapta a la pantalla gràcies a max-width: 100%)
    'style_nav_header_background': '#0077C8',
    'prev_next_buttons_location': 'bottom',
    'titles_only': False,
}

# Metadades HTML
html_title = f'{project} v{release}'
html_short_title = project
html_show_sourcelink = False
html_show_sphinx = False

# CSS personalitzat per fer el contingut responsiu (ocupa tota l'amplada disponible)
html_static_path = ['_static']
html_css_files = ['custom.css']

# ---------------------------------------------------------------------------
# Formats de font acceptats
# ---------------------------------------------------------------------------
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
} 