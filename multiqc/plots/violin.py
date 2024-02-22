""" MultiQC functions to plot a violin plot """

import logging
from typing import List, Dict, Optional, Union

from multiqc.plots import table_object
from multiqc.utils import config
from multiqc.plots.plotly import violin

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


def plot(data: List[Dict], headers: Optional[Union[List[Dict], Dict]] = None, pconfig=None):
    """Helper HTML for a violin plot.
    :param data: A list of data dicts
    :param headers: A list of dicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :return: HTML string
    """
    # Make a datatable object
    dt = table_object.DataTable(data, headers, pconfig)

    mod = get_template_mod()
    if "violin" in mod.__dict__ and callable(mod.violin):
        # noinspection PyBroadException
        try:
            return mod.violin(dt)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return violin.plot(dt)
