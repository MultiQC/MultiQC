"""
========
 plotly
========

The Plotly rendering template: the former 'default' template, kept for reports
that want interactive Plotly plots. Requires the optional `plotly` dependency
(`pip install 'multiqc[plotly]'`). Selected via `multiqc --template plotly`.

The main 'default' template now renders plots with Apache ECharts, which needs no
extra dependencies for interactive or static export.
"""

import os

template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Template configuration - overrides user config
plotting_engine = "plotly"
