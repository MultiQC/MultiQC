"""
========
 simple
========

This theme attempts to generate a report which is a simple
as possible - no JavaScript where possible. The resulting
report is hopefully suitable for e-mailing, converting to PDF
and printing.

"""

import os
import importlib

from multiqc import config

template_parent = "original"

config.plots_force_flat = True
config.simple_output = True

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Import template functions from parent
try:
    parent_mod = importlib.import_module(f"multiqc.templates.{template_parent}")
    if hasattr(parent_mod, "template_functions"):
        template_functions = parent_mod.template_functions
except ImportError:
    pass
