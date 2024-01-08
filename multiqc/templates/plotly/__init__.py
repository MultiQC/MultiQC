"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os
from .plots import bar, line, scatter, heatmap, violin
from ...utils import config

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

bargraph = bar.plot
linegraph = line.plot
scatter = scatter.plot
heatmap = heatmap.plot
beeswarm = violin.plot
violin = violin.plot

__all__ = ["bargraph", "linegraph", "scatter", "heatmap", "beeswarm", "violin"]

if config.development:
    copy_files = ["assets"]
