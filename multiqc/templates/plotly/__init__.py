"""
========
 plotly
========

The Plotly rendering template, kept for reports that want interactive Plotly plots.
Requires the optional `plotly` dependency (`pip install 'multiqc[plotly]'`). Selected
via `multiqc --template plotly`.

A thin child of the `default` template: it inherits every shared HTML file (base, nav,
toolbox, head, includes, footer, ...) from default and overrides only the plotting layer,
its own `src/js`/`compiled` Plotly renderer bundle and the vendored Plotly asset. The
shared `includes.html`/`footer.html` switch the plot bundle and credit on
`config.plotting_engine`, which this template sets to "plotly".
"""

import os

template_dir = os.path.dirname(__file__)
template_parent = "default"
base_fn = "base.html"

# Template configuration - overrides user config
plotting_engine = "plotly"
