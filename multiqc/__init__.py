#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from collections import OrderedDict
import json
import os
import random

letters = 'abcdefghijklmnopqrstuvwxyz'

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
        s_name = s_name.split("_tophat",1)[0]
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
        if config.get('id') is None:
            config['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
        
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
    
    def plot_bargraph (self, data, cats=None, config={}):
        """ Plot a horizontal bar graph. Expects a 2D dict of sample
        data. Also can take info about categories. There are quite a
        few variants of how to use this function, see CONTRIBUTING.md
        for documentation and examples. """
        
        # Not given any cats - find them from the data
        if cats is None:
            cats = list(set(k for s in data.keys() for k in data[s].keys() ))
        
        # Given a list of cats - turn it into a dict
        if type(cats) is list:
            newcats = OrderedDict()
            for c in cats:
                newcats[c] = {'name': c}
            cats = newcats
        
        # Parse the data into a HighCharts friendly format
        hc_samples = data.keys()
        hc_data = list()
        for c in cats.keys():
            thisdata = list()
            for s in data.keys():
                thisdata.append(data[s][c])
            if max(thisdata) > 0:
                thisdict = { 'name': cats[c]['name'], 'data': thisdata }
                if 'color' in cats[c]:
                    thisdict['color'] = cats[c]['color']
                hc_data.append(thisdict)
        
        # Build the HTML
        if config.get('id') is None:
            config['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
        html = ''
        
        # Counts / Percentages Switch
        if config.get('cpswitch') is not False:
            c_active = 'active' if config.get('cpswitch_c_active', True) is True else ''
            p_active = 'active' if c_active is not 'active' else ''
            c_label = config.get('cpswitch_counts_label', 'Counts')
            p_label = config.get('cpswitch_percent_label', 'Percentages')
            html += '<div class="btn-group switch_group"> \n\
    			<button class="btn btn-default btn-sm {c_a}" data-action="set_numbers" data-target="#{id}">{c_l}</button> \n\
    			<button class="btn btn-default btn-sm {p_a}" data-action="set_percent" data-target="#{id}">{p_l}</button> \n\
    		</div>'.format(id=config['id'], c_a=c_active, p_a=p_active, c_l=c_label, p_l=p_label)
        
        # Plot and javascript function
        html += '<div id="{id}" class="hc-plot"></div> \n\
        <script type="text/javascript"> \n\
            $(function () {{ plot_stacked_bar_graph("#{id}", {s}, {d}, {c}); }}); \
        </script>'.format(id=config['id'], s=json.dumps(hc_samples), d=json.dumps(hc_data), c=json.dumps(config));
        
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
