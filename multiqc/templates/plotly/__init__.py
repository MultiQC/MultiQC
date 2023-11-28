"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os

from .plots import bar_plot, line_plot, scatter_plot

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

base64_plots = False


bargraph = bar_plot.bar_plot
linegraph = line_plot.line_plot
scatter_plot = scatter_plot.scatter_plot

__all__ = ["bargraph", "linegraph", "scatter_plot"]
