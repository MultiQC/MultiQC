#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from collections import OrderedDict
import json
import os
import random

from multiqc import config

letters = 'abcdefghijklmnopqrstuvwxyz'

class BaseMultiqcModule(object):

    def __init__(self):
        pass

    def clean_s_name(self, s_name, root, trimmed=True):
        """ Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary.
        :param s_name: The sample name to clean
        :param root: The directory path that this file is within
        :param prepend_dirs: boolean, whether to prepend dir name to s_name
        :param trimmed: boolean, remove common trimming suffixes from name?
        :return: The cleaned sample name, ready to be used
        """
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
        if config.prepend_dirs:
            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
        return s_name
    
    
    def plot_xy_data(self, data, config={}, original_plots=[]):
        """ Plot a line graph with X,Y data. See CONTRIBUTING.md for
        further instructions on use.
        :param data: 2D dict, first keys as sample names, then x:y data pairs
        :param original_plots: optional list of dicts with keys 's_name' and 'img_path'
        :param config: optional dict with config key:value pairs. See CONTRIBUTING.md
        :param original_plots: optional list specifying original plot images. Each dict
                               should have a key 's_name' and 'img_path'
        :return: HTML and JS, ready to be inserted into the page
        """
        
        # Given one dataset - turn it into a list
        if type(data) is not list:
            data = [data]
        
        # Generate the data dict structure expected by HighCharts series
        plotdata = list()
        for d in data:
            thisplotdata = list()
            for s in sorted(d.keys()):
                pairs = list()
                for k, p in d[s].items():
                    pairs.append([k, p])
                thisplotdata.append({
                    'name': s,
                    'data': pairs
                })
            plotdata.append(thisplotdata)
        
        # Build the HTML for the page
        if config.get('id') is None:
            config['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
        html = ''
        
        # Buttons to cycle through different datasets
        if len(plotdata) > 1:
            html += '<div class="btn-group switch_group">\n'
            for k, p in enumerate(plotdata):
                active = 'active' if k == 0 else ''
                try: name = config['data_labels'][k]['name']
                except: name = k+1
                try: ylab = 'data-ylab="{}"'.format(config['data_labels'][k]['ylab'])
                except: ylab = 'data-ylab="{}"'.format(name) if name != k+1 else ''
                html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} data-newdata="{id}_datasets[{k}]" data-target="#{id}">{n}</button>\n'.format(a=active, id=config['id'], n=name, y=ylab, k=k)
            html += '</div>\n\n'
        
        # Markup needed if we have the option of clicking through to original plot images
        if len(original_plots) > 0:
            config['tt_label'] = 'Click to show original plot.<br>{}'.format(config.get('tt_label', '{point.x}'))
            if len(original_plots) > 1:
                next_prev_buttons = '<div class="clearfix"><div class="btn-group btn-group-sm"> \n\
                    <a href="#{prev}" class="btn btn-default original_plot_prev_btn" data-target="#{id}">&laquo; Previous</a> \n\
                    <a href="#{next}" class="btn btn-default original_plot_nxt_btn" data-target="#{id}">Next &raquo;</a> \n\
                </div></div>'.format(id=config['id'], prev=original_plots[-1]['s_name'], next=original_plots[1]['s_name'])
            else:
                next_prev_buttons = ''
            html += '<p class="text-muted instr">Click to show original FastQC plot.</p>\n\
                    <div id="fastqc_quals" class="hc-plot-wrapper"> \n\
                        <div class="showhide_orig" style="display:none;"> \n\
                            <h4><span class="s_name">{n}</span></h4> \n\
                            {b} <img data-toggle="tooltip" title="Click to return to overlay plot" class="original-plot" src="{f}"> \n\
                        </div>\n\
                        <div id="{id}" class="hc-plot"></div> \n\
                    </div>'.format(id=config['id'], b=next_prev_buttons, n=original_plots[0]['s_name'], f=original_plots[0]['img_path'])
            orig_plots = 'var {id}_orig_plots = {d}; \n'.format(id=config['id'], d=json.dumps(original_plots))
            config['orig_click_func'] = True # Javascript prints the click function
        
        # Regular plots (no original images)
        else:
            html += '<div id="{id}" class="hc-plot"></div> \n'.format(id=config['id'])
            orig_plots = ''
        
        # Javascript with data dump
        html += '<script type="text/javascript"> \n\
            var {id}_datasets = {d}; \n\
            {o} \
            $(function () {{ plot_xy_line_graph("#{id}", {id}_datasets[0], {c}); }}); \n\
        </script>'.format(id=config['id'], d=json.dumps(plotdata), c=json.dumps(config), o=orig_plots);
        return html
    
    
    def plot_bargraph (self, data, cats=None, config={}):
        """ Plot a horizontal bar graph. Expects a 2D dict of sample
        data. Also can take info about categories. There are quite a
        few variants of how to use this function, see CONTRIBUTING.md
        for documentation and examples.
        :param data: 2D dict, first keys as sample names, then x:y data pairs
        :param cats: optnal list, dict or OrderedDict with plot categories
        :param config: optional dict with config key:value pairs
        :return: HTML and JS, ready to be inserted into the page
        """
        
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
        hc_samples = sorted(list(data.keys()))
        hc_data = list()
        for c in cats.keys():
            thisdata = list()
            for s in hc_samples:
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
            if config.get('cpswitch_c_active', True) is True:
                c_active = 'active'
                p_active = ''
            else:
                c_active = ''
                p_active = 'active'
                config['stacking'] = 'percent'
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
        """ Converts a dict to a CSV string
        :param d: 2D dictionary, first keys sample names and second key
                  column headers
        :param delim: optional delimiter character. Default: \t
        :return: Flattened string, suitable to write to a CSV file.
        """

        h = None # We make a list of keys to ensure consistent order
        l = list()
        for sn in d:
            if h is None:
                h = list(d[sn].keys())
                l.append(delim.join([''] + h))
            thesefields = [sn] + [ str(d[sn].get(k, '')) for k in h ]
            l.append( delim.join( thesefields ) )
        return ('\n'.join(l)).encode('utf-8', 'ignore').decode('utf-8')
