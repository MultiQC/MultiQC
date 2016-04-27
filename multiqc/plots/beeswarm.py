#!/usr/bin/env python

""" MultiQC functions to plot a beeswarm group """

import json
import logging
import os
import random

from multiqc.utils import config
from multiqc.plots import table_object

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=[], pconfig={}):
    """ Helper HTML for a beeswarm plot.
    :param data: A list of data dicts
    :param headers: A list of Dicts / OrderedDicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :return: HTML string
    """
    
    # Make a datatable object
    dt = table_object.datatable(data, headers, pconfig)
    
    return make_beeswarm_plot( dt )
    
    
def make_beeswarm_plot(dt):    
    categories = []
    s_names = []
    data = []
    for idx, d in dt.data.items():
    
        for k in headers[idx].keys():
            
            rid = headers[idx][k]['rid']
            bcol = 'rgb({})'.format(headers[idx][k].get('colour', '204,204,204'))

            categories.append({
                'title': headers[idx][k]['title'],
                'description': headers[idx][k]['description'],
                'max': headers[idx][k]['dmax'],
                'min': headers[idx][k]['dmin'],
                'suffix': headers[idx][k].get('suffix', ''),
                'decimalPlaces': headers[idx][k].get('decimalPlaces', '2'),
                'bordercol': bcol
            });
            
            # Add the data
            thisdata = []
            these_snames = []
            for (sname, samp) in d:
                if k in samp:
                    
                    val = samp[k]
                    general_stats_raw[sname][rid] = val
                    
                    if 'modify' in headers[idx][k] and callable(headers[idx][k]['modify']):
                        val = headers[idx][k]['modify'](val)
                    
                    thisdata.append(val)
                    these_snames.append(sname)
            
            data.append(thisdata)
            s_names.append(these_snames)
    
    # Plot and javascript function
    global general_stats_beeswarm_html
    general_stats_beeswarm_html = '<div class="hc-plot-wrapper"><div id="general_stats_beeswarm" class="hc-plot not_rendered hc-beeswarm-plot"><small>loading..</small></div></div> \n\
    <script type="text/javascript"> \n\
        mqc_plots["general_stats_beeswarm"] = {{ \n\
            "plot_type": "beeswarm", \n\
            "samples": {s}, \n\
            "datasets": {d}, \n\
            "categories": {c} \n\
        }} \n\
    </script>'.format(s=json.dumps(s_names), d=json.dumps(data), c=json.dumps(categories));