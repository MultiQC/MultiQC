#!/usr/bin/env python

""" MultiQC functions to plot a scatter plot """

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


def plot(plotdata, pconfig=None):
    return highcharts_scatter_plot(plotdata, pconfig)


def highcharts_scatter_plot(plotdata, pconfig=None):
    """
    Build the HTML needed for a HighCharts scatter plot. Should be
    called by scatter.plot(), which properly formats input data.
    """
    if pconfig is None:
        pconfig = {}

    # Get the plot ID
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_hcplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    # Build the HTML for the page
    html = '<div class="mqc_hcplot_plotgroup">'

    # Buttons to cycle through different datasets
    if len(plotdata) > 1:
        html += '<div class="btn-group hc_switch_group">\n'
        for k, p in enumerate(plotdata):
            active = "active" if k == 0 else ""
            dls = [dl if isinstance(dl, dict) else {"name": dl} for dl in pconfig.get("data_labels", [])]
            try:
                name = dls[k]["name"]
            except Exception:
                name = k + 1
            try:
                ylab = f"data-ylab=\"{dls[k]['ylab']}\""
            except Exception:
                ylab = f'data-ylab="{name}"' if name != k + 1 else ""
            try:
                ymax = f"data-ymax=\"{dls[k]['ymax']}\""
            except Exception:
                ymax = ""
            try:
                xlab = f"data-xlab=\"{dls[k]['xlab']}\""
            except Exception:
                xlab = f'data-xlab="{name}"' if name != k + 1 else ""
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} {xl} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(
                a=active, id=pconfig["id"], n=name, y=ylab, ym=ymax, xl=xlab, k=k
            )
        html += "</div>\n\n"

    # The plot div
    html += '<div class="hc-plot-wrapper"{height}><div id="{id}" class="hc-plot not_rendered hc-scatter-plot"><small>loading..</small></div></div></div> \n'.format(
        id=pconfig["id"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    # Reverse order of dots in plotdata as the z-order is reversed in Highcharts compared to Plotly
    for d in plotdata:
        d.reverse()

    report.plot_data[pconfig["id"]] = {"plot_type": "scatter", "datasets": plotdata, "config": pconfig}

    return html
