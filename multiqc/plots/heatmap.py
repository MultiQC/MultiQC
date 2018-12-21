#!/usr/bin/env python

""" MultiQC functions to plot a heatmap """

from __future__ import print_function
import logging
import random

from multiqc.utils import config, report

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, xcats, ycats=None, pconfig=None):
    """ Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values.
    :param xcats: Labels for x axis
    :param ycats: Labels for y axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """

    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if 'id' in pconfig and pconfig['id'] and pconfig['id'] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig['id']].items():
            pconfig[k] = v

    if ycats is None:
        ycats = xcats

    # Make a plot
    return highcharts_heatmap(data, xcats, ycats, pconfig)



def highcharts_heatmap (data, xcats, ycats, pconfig=None):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by plot_xy_data, which properly formats input data.
    """
    if pconfig is None:
        pconfig = {}

    # Reformat the data for highcharts
    pdata = []
    for i, arr in enumerate(data):
        for j, val in enumerate(arr):
            pdata.append([j,i,val])

    # Get the plot ID
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig['id'] = report.save_htmlid(pconfig['id'])

    # Build the HTML for the page
    html = '<div class="mqc_hcplot_plotgroup">'

    # The 'sort by highlights button'
    html += """<div class="btn-group hc_switch_group">
        <button type="button" class="mqc_heatmap_sortHighlight btn btn-default btn-sm" data-target="#{id}" disabled="disabled">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
    </div>""".format(id=pconfig['id'])

    # The plot div
    html += '<div class="hc-plot-wrapper"><div id="{id}" class="hc-plot not_rendered hc-heatmap"><small>loading..</small></div></div></div> \n'.format(id=pconfig['id'])

    report.num_hc_plots += 1

    report.plot_data[pconfig['id']] = {
        'plot_type': 'heatmap',
        'data': pdata,
        'xcats': xcats,
        'ycats': ycats,
        'config': pconfig
    }

    return html
