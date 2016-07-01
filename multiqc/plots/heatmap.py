#!/usr/bin/env python

""" MultiQC functions to plot a heatmap """

from __future__ import print_function
import json
import logging
import os
import random

from multiqc.utils import config

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

# Load the template so that we can access it's configuration
template_mod = config.avail_templates[config.template].load()

def plot (data, xcats, ycats=None, pconfig={}):
    """ Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values.
    :param xcats: Labels for x axis
    :param ycats: Labels for y axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """
    
    if ycats is None:
        ycats = xcats
    
    # Make a plot
    return highcharts_heatmap(data, xcats, ycats, pconfig)



def highcharts_heatmap (data, xcats, ycats, pconfig={}):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by plot_xy_data, which properly formats input data.
    """
    
    # Reformat the data for highcharts
    pdata = []
    for i, arr in enumerate(data):
        for j, val in enumerate(arr):
            pdata.append([j,i,val])
    
    # Build the HTML for the page
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
    html = '<div class="mqc_hcplot_plotgroup">'
    
    # The plot div
    html += '<div class="hc-plot-wrapper"><div id="{id}" class="hc-plot not_rendered hc-heatmap"><small>loading..</small></div></div></div> \n'.format(id=pconfig['id'])
    
    # Javascript with data dump
    html += '<script type="text/javascript"> \n\
        mqc_plots["{id}"] = {{ \n\
            "plot_type": "heatmap", \n\
            "data": {d}, \n\
            "xcats": {x}, \n\
            "ycats": {y}, \n\
            "config": {c} \n\
        }} \n\
    </script>'.format(id=pconfig['id'], d=json.dumps(pdata), x=json.dumps(xcats), y=json.dumps(ycats), c=json.dumps(pconfig));
    
    return html

