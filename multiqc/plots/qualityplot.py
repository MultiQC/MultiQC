#!/usr/bin/env python

""" Functions to plot miniature boxplots for ber pase quality """

import sys
import logging


logger = logging.getLogger(__name__)





try:
    # Import matplot lib but avoid default X environment
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logger.debug("Using matplotlib version {}".format(matplotlib.__version__))
except Exception as e:
    # MatPlotLib can break in a variety of ways. Fake an error message and continue without it if so.
    # The lack of the library will be handled when plots are attempted
    print("##### ERROR! MatPlotLib library could not be loaded!    #####", file=sys.stderr)
    print("##### Flat plots will instead be plotted as interactive #####", file=sys.stderr)
    print(e)






letters = "abcdefghijklmnopqrstuvwxyz"
# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None
def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod




def plot(data, pconfig=None):
    """Plot per base quality as miniature boxplots.
    :param data: dictionary. Keys are samples. Values are dictionaries of quality statistics per base.
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """



