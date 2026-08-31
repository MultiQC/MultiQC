"""
========
 echarts
========

Back-compat alias for `multiqc --template echarts`. The main `default` template now
renders with Apache ECharts, so this is a thin child of `default` that inherits all of
its files (includes.html, compiled JS/CSS, the vendored ECharts bundle). Kept so existing
`--template echarts` invocations keep working; new reports should just use `default`.
"""

import os

template_parent = "default"
template_dir = os.path.dirname(__file__)
base_fn = "base.html"

# Inherited from `default` anyway; set explicitly so the engine resolves even if the
# parent chain changes.
plotting_engine = "echarts"
