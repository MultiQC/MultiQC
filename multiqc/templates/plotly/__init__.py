"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os

from .plots import bar, line, scatter

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

base64_plots = False


bargraph = bar.plot
linegraph = line.plot
scatter = scatter.plot

__all__ = ["bargraph", "linegraph", "scatter"]
