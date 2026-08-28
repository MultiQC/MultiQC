"""
========
 echarts
========

A child theme of 'default' that renders plots with Apache ECharts instead of
Plotly. Selected via `multiqc --template echarts`.

"""

import importlib
import os

template_parent = "default"
template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Template configuration - overrides user config
plotting_engine = "echarts"

# Import template functions from parent
try:
    parent_mod = importlib.import_module(f"multiqc.templates.{template_parent}")
    if hasattr(parent_mod, "template_functions"):
        template_functions = parent_mod.template_functions
except ImportError:
    pass
