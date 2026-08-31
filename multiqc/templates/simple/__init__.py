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

from multiqc import config

template_parent = "original"

# Render flat/static images via the ECharts SSR engine (browserless): the simple
# template emits no interactive JS, so its plots are always static images, and ECharts
# owns all static image export. This also fixes `--pdf` (which forces the simple template).
plotting_engine = "echarts"

config.plots_force_flat = True
config.simple_output = True

template_dir = os.path.dirname(__file__)
base_fn = "base.html"
