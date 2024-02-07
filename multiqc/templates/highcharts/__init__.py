"""
=============
 highcharts
=============

Legacy template using HighCharts and matplotlib as a plotting backend.

"""
import os
from .plots import bargraph, beeswarm, linegraph, scatter, heatmap, table

from multiqc.utils import config

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

bargraph = bargraph.plot
linegraph = linegraph.plot
scatter = scatter.plot
heatmap = heatmap.plot
violin = beeswarm.plot
beeswarm = beeswarm.plot
table = table.plot

__all__ = ["bargraph", "linegraph", "scatter", "heatmap", "beeswarm", "violin", "table"]

if config.development:
    copy_files = ["assets"]
