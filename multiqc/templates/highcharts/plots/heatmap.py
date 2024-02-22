#!/usr/bin/env python

""" MultiQC functions to plot a heatmap """


import logging
import random

from multiqc.utils import config, report

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


def plot(data, xcats, ycats, pconfig=None):
    """Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values.
    :param xcats: Labels for x axis
    :param ycats: Labels for y axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """
    return highcharts_heatmap(data, xcats, ycats, pconfig)


def highcharts_heatmap(data, xcats, ycats, pconfig=None):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by plot_xy_data, which properly formats input data.
    """
    if pconfig is None:
        pconfig = {}

    # Reformat the data for highcharts
    pdata = []
    minval = None
    maxval = None
    for i, arr in enumerate(data):
        for j, val in enumerate(arr):
            pdata.append([j, i, val])
            if isinstance(val, (int, float)):
                if minval is None or val < minval:
                    minval = val
                if maxval is None or val > maxval:
                    maxval = val

    if "min" not in pconfig:
        pconfig["min"] = minval
    if "max" not in pconfig:
        pconfig["max"] = maxval

    # Get the plot ID
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_hcplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    # Build the HTML for the page
    html = """
    <div class="mqc_hcplot_plotgroup">
        <div class="mqc_hcplot_range_sliders">
            <div>
                <label for="{id}_range_slider_min_txt">Min:</label>
                <input id="{id}_range_slider_min_txt" type="number" class="form-control" value="{min}" min="{min}" max="{max}" data-minmax="min" data-target="{id}" />
                <input id="{id}_range_slider_min" type="range" value="{min}" min="{min}" max="{max}" step="any" data-minmax="min" data-target="{id}" />
            </div>
            <div>
                <label for="{id}_range_slider_max_txt">Max:</label>
                <input id="{id}_range_slider_max_txt" type="number" class="form-control" value="{max}" min="{min}" max="{max}" data-minmax="max" data-target="{id}" />
                <input id="{id}_range_slider_max" type="range" value="{max}" min="{min}" max="{max}" step="any" data-minmax="max" data-target="{id}" />
            </div>
        </div>
        <div class="hc-plot-wrapper"{height}>
            <div id="{id}" class="hc-plot not_rendered hc-heatmap">
                <small>loading..</small>
            </div>
        </div>
    </div> \n""".format(
        id=pconfig["id"],
        min=pconfig["min"],
        max=pconfig["max"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[pconfig["id"]] = {
        "plot_type": "heatmap",
        "data": pdata,
        "xcats": xcats,
        "ycats": ycats,
        "config": pconfig,
    }

    return html
