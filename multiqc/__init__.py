#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

import json
import os
import random
import string

class BaseMultiqcModule(object):

    def __init__(self):
        pass

    def clean_s_name(self, s_name, root, prepend_dirs=False, trimmed=True):
        """ Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary. """
        # Split then take first section to remove everything after these matches
        s_name = s_name.split(".gz",1)[0]
        s_name = s_name.split(".fastq",1)[0]
        s_name = s_name.split(".fq",1)[0]
        s_name = s_name.split(".bam",1)[0]
        s_name = s_name.split(".sam",1)[0]
        if trimmed:
            s_name = s_name.split("_val_1",1)[0]
            s_name = s_name.split("_val_2",1)[0]
            s_name = s_name.split("_trimmed",1)[0]
        if prepend_dirs:
            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
        return s_name
    
    def plot_xy_data(self, data, config={}):
        """ Plot a line graph with X,Y data. Expects a 2D dict with 
        first keys as sample names, then x: y data pairs.
        See CONTRIBUTING.md for further instructions on use. """
        
        # Generate a random HTML id if not given
        if config['id'] is None:
            config['id'] = ''.join(random.sample(string.lowercase+string.digits, 10))
        
        # Generate the data dict structure expected by HighCharts series
        plotdata = list()
        for s in sorted(data):
            pairs = list()
            for k, p in iter(sorted(data[s].items())):
                pairs.append([k, p])
            plotdata.append({
                'name': s,
                'data': pairs
            })
        
        # Build the HTML for the page
        html = '<div id="{n}" class="hc-plot"></div> \n\
        <script type="text/javascript"> \n\
            $(function () {{ plot_xy_line_graph("#{n}", {d}, {c}); }}); \n\
        </script>'.format(n=config['id'], d=json.dumps(data), c=json.dumps(config));
        return html
    
    def dict_to_csv (self, d, delim="\t"):
        """ Takes a 2D dictionary and returns a string suitable for
        writing to a .csv file. First key should be sample name
        (row header), second key should be field (column header)
        Takes dictionary as input, optional 'delim' field to specify
        column delimiter (default: tab) """

        h = None # We make a list of keys to ensure consistent order
        l = list()
        for sn in d:
            if h is None:
                h = list(d[sn].keys())
                l.append(delim.join([''] + h))
            thesefields = [sn] + [ str(d[sn].get(k, '')) for k in h ]
            l.append( delim.join( thesefields ) )
        return ('\n'.join(l)).encode('utf-8', 'ignore').decode('utf-8')
