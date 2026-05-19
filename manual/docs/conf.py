import sys
import os
sys.path.insert(0, os.path.abspath('.'))

project = 'Manuals'
author = 'Carlos Herrera Vázquez & Alexandre Gràcia-Calvo & Eduardo Prieto-Araujo'
release = '0.1'

extensions = [
    'myst_parser',
]

html_theme = 'sphinx_rtd_theme'

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
} 