"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os
from .plots import bar, line, scatter, heatmap, violin, table
from ...utils import config

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

bargraph = bar.plot
linegraph = line.plot
scatter = scatter.plot
heatmap = heatmap.plot
violin = violin.plot
table = table.plot

__all__ = ["bargraph", "linegraph", "scatter", "heatmap", "violin", "table"]

if config.development:
    copy_files = ["assets"]
