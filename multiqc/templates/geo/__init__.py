"""
=====
 Geo
=====

Welcome to Geo - an example template / easter egg.

Theme courtesy of the wonderful Geo for Bootstrap theme:
http://code.divshot.com/geo-bootstrap/

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
