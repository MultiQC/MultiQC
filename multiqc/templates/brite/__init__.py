"""
========
 brite
========

A child theme of 'default' using the Bootswatch Brite theme.

Brite is a bright, high-contrast theme from Bootswatch.

"""

import os
import importlib

template_parent = "default"
template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Import template functions from parent
try:
    parent_mod = importlib.import_module(f"multiqc.templates.{template_parent}")
    if hasattr(parent_mod, "template_functions"):
        template_functions = parent_mod.template_functions
except ImportError:
    pass
