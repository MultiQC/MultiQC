#!/usr/bin/env python

""" MultiQC functions to plot a scatter plot """

import json
import logging
import os
import random

from multiqc.utils import config

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

# Load the template so that we can access it's configuration
template_mod = config.avail_templates[config.template].load()

def plot (data, pconfig={}):
    """ Plot a scatter plot with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    
    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]
    
    # Generate the data dict structure expected by HighCharts series
    plotdata = list()
    for ds in data:
        d = list()
        for s_name in ds:
            for k in ds[s_name]:
                this_series = { 'name': s_name, 'x': k['x'], 'y': k['y'] }
                try:
                    this_series['color'] = pconfig['colors'][s_name]
                except: pass
                d.append(this_series)
        plotdata.append(d)
    
    # Add on annotation data series
    try:
        for s in pconfig['extra_series']:
            plotdata[0].append(s)
    except KeyError:
        pass
    
    # Make a plot
    return highcharts_scatter_plot(plotdata, pconfig)

def highcharts_scatter_plot (plotdata, pconfig={}):
    """
    Build the HTML needed for a HighCharts scatter plot. Should be
    called by scatter.plot(), which properly formats input data.
    """
    
    # Build the HTML for the page
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
    html = '<div class="mqc_hcplot_plotgroup">'
    
    # Buttons to cycle through different datasets
    if len(plotdata) > 1:
        html += '<div class="btn-group hc_switch_group">\n'
        for k, p in enumerate(plotdata):
            active = 'active' if k == 0 else ''
            try: name = pconfig['data_labels'][k]['name']
            except: name = k+1
            try: ylab = 'data-ylab="{}"'.format(pconfig['data_labels'][k]['ylab'])
            except: ylab = 'data-ylab="{}"'.format(name) if name != k+1 else ''
            try: ymax = 'data-ymax="{}"'.format(pconfig['data_labels'][k]['ymax'])
            except: ymax = ''
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(a=active, id=pconfig['id'], n=name, y=ylab, ym=ymax, k=k)
        html += '</div>\n\n'
    
    # The plot div
    html += '<div class="hc-plot-wrapper"><div id="{id}" class="hc-plot not_rendered hc-scatter-plot"><small>loading..</small></div></div></div> \n'.format(id=pconfig['id'])
    
    # Javascript with data dump
    html += '<script type="text/javascript"> \n\
        mqc_plots["{id}"] = {{ \n\
            "plot_type": "scatter", \n\
            "datasets": {d}, \n\
            "config": {c} \n\
        }} \n\
    </script>'.format(id=pconfig['id'], d=json.dumps(plotdata), c=json.dumps(pconfig));
    
    return html
    