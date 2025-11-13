"""
========
 crazy
========

A child theme of 'default' - experimental / testing template.

"""

import os
import importlib

template_parent = "default"
template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Template configuration - overrides user config
template_dark_mode = False  # Disable dark mode toggle for this template
plot_font_family = '"Rock Salt", cursive'

# Import template functions from parent
try:
    parent_mod = importlib.import_module(f"multiqc.templates.{template_parent}")
    if hasattr(parent_mod, "template_functions"):
        template_functions = parent_mod.template_functions
except ImportError:
    pass
