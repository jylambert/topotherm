# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent  # docs/source
DOCS_DIR = THIS_DIR.parent  # docs
ROOT = DOCS_DIR.parent  # repo root

# Add the root directory of your project to sys.path
sys.path.insert(0, str(ROOT))


project = "topotherm"
copyright = "2025, Lambert, Jerry and Ceruti, Amedeo"
author = "Jerry Lambert & Amedeo Ceruti"
release = "0.6.0"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'autoapi.extension',  # Must be included
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",  # automatic summary tables
    "sphinx.ext.napoleon",  # numpy/google style
    "sphinx_autodoc_typehints",  # show type hints in docs
    "sphinx.ext.viewcode",  # add "view source"
    "sphinx.ext.intersphinx",  # link to external docs
    "sphinx_copybutton",  # copy button in code blocks
    # 'myst_parser', we need this for markdown, if all works delete
    "autoapi.extension",
    "sphinx.ext.mathjax",
]
# Prefer NumPy-style sections
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True  # show "Parameters"
napoleon_use_rtype = True  # show "Returns"
napoleon_attr_annotations = True
# Nice rendering tweaks
autodoc_typehints = "description"  # move type hints under params/returns
always_document_param_types = True  # Sphinx ≥ 7

autosummary_generate = True
# myst_enable_extensions = [
#     "colon_fence",
#     "deflist",
#     "linkify",
# ]

# make clear we only want to use .rst and index
source_suffix = ".rst"
root_doc = "index"

templates_path = ["_templates"]
exclude_patterns = []

# AutoAPI configuration (this is what keeps your function/class list fresh)
# -----------------------------------------------------------------------------
autoapi_type = "python"
# Point directly at your package in the repo root:
autoapi_dirs = [str(ROOT / "topotherm")]
autoapi_root = "autoapi"
autoapi_add_toctree_entry = True
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]

# -----------------------------------------------------------------------------
# HTML output
# -----------------------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']

html_theme_options = {
    "collapse_navigation": False,  # don't collapse siblings when navigating
    "sticky_navigation": True,  # keep sidebar pinned while scrolling
    "navigation_depth": 5,  # how deep to expand
    "titles_only": False,  # show sections under pages
}
# -----------------------------------------------------------------------------
# Intersphinx (optional; keep if you plan to cross-link to these)
# -----------------------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}
