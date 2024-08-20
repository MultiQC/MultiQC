"""MultiQC functions to plot a heatmap"""

import logging
from typing import Union, Dict, Optional, List

from multiqc import config
from multiqc.plots.plotly import heatmap
from multiqc.plots.plotly.heatmap import HeatmapConfig

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(
    data,
    xcats: Optional[List[Union[str, int]]] = None,
    ycats: Optional[List[Union[str, int]]] = None,
    pconfig: Union[Dict, HeatmapConfig, None] = None,
) -> Union[heatmap.HeatmapPlot, str]:
    """Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values; or a dict of dicts
    :param xcats: Labels for x-axis
    :param ycats: Labels for y-axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """
    pconf: HeatmapConfig = HeatmapConfig.from_pconfig_dict(pconfig)

    if ycats is None:
        ycats = xcats

    # Make a plot
    mod = get_template_mod()
    if "heatmap" in mod.__dict__ and callable(mod.heatmap):
        # noinspection PyBroadException
        try:
            return mod.heatmap(data, xcats, ycats, pconf)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return heatmap.plot(data, pconf, xcats, ycats)
