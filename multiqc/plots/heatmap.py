""" MultiQC functions to plot a heatmap """


import logging

from multiqc.utils import config
from multiqc.plots.plotly import heatmap

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(data, xcats=None, ycats=None, pconfig=None):
    """Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values; or a dict of dicts
    :param xcats: Labels for x-axis
    :param ycats: Labels for y-axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

    if ycats is None:
        ycats = xcats

    # Make a plot
    mod = get_template_mod()
    if "heatmap" in mod.__dict__ and callable(mod.heatmap):
        try:
            return mod.heatmap(data, xcats, ycats, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return heatmap.plot(data, pconfig, xcats, ycats)
