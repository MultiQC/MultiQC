"""
=========
 gathered
=========

This theme extends 'default'.  Unlike 'default', it organizes all sections
under the 'Detailed Statistics' header, rather than grouping them by
module.
"""

import os
import importlib

template_dir = os.path.dirname(__file__)
template_parent = "default"
base_fn = "base.html"

# Import template functions from parent
try:
    parent_mod = importlib.import_module(f"multiqc.templates.{template_parent}")
    if hasattr(parent_mod, "template_functions"):
        template_functions = parent_mod.template_functions
except ImportError:
    pass
