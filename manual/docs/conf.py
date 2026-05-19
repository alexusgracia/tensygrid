import sys
import os
sys.path.insert(0, os.path.abspath('.'))

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