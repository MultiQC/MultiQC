"""
=============
 plolty
=============

Theme to work on porting the default plotting library
from HighCharts to Plotly.

"""
import os

from .plots import bargraph

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"


def custom_linegraph(plot_data, pconfig):
    return "<h1>Awesome line graph here</h1>"


# linegraph = custom_linegraph
bargraph = bargraph.plotly_bargraph
