"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os

from .plots import bargraph, linegraph

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

base64_plots = False


bargraph = bargraph.bargraph
linegraph = linegraph.linegraph
