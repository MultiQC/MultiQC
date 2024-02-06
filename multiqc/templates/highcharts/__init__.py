"""
=============
 highcharts
=============

Legacy template using HighCharts and matplotlib as a plotting backend.

"""
import os

from multiqc.utils import config

template_parent = "default"

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

if config.development:
    copy_files = ["assets"]
