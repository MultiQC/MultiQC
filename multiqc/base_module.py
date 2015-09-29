#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

from __future__ import print_function
from collections import OrderedDict
import io
import json
import mimetypes
import os
import random
import logging

from multiqc import config
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

class BaseMultiqcModule(object):

    def __init__(self, name='base', anchor='base', target='',href='', info='', extra=''):
        self.name = name
        self.anchor = anchor
        if not target:
            target = self.name
        self.intro = '<p><a href="{0}" target="_blank">{1}</a> {2}</p>{3}'.format(
            href, target, info, extra
        )

    def find_log_files(self, fn_match=None, contents_match=None, filehandles=False):
        """
        Search the analysis directory for log files of interest. Can take either a filename
        suffix or a search string to return only log files that contain relevant info.
        :param fn_match: Optional string or list of strings. Filename suffixes to search for.
        :param contents_match: Optional string or list of strings to look for in file.
        NB: Both searches return file if *any* of the supplied strings are matched.
        :param filehandles: Set to true to return a file handle instead of slurped file contents
        :return: Yields a set with two items - a sample name generated from the filename
                 and either the file contents or file handle for the current matched file.
                 As yield is used, the function can be iterated over without 
        """
        for directory in config.analysis_dir:
            for root, dirnames, filenames in os.walk(directory, followlinks=True):
                
                for fn in filenames:
                    
                    # Ignore files set in config
                    if fn in config.fn_ignore_files:
                        continue
                    
                    # Make a sample name from the filename
                    s_name = self.clean_s_name(fn, root)
                    
                    # Make search strings into lists if a string is given
                    if type(fn_match) is str:
                        fn_match = [fn_match]
                    if type(contents_match) is str:
                        contents_match = [contents_match]
                    
                    # Search for file names ending in a certain string
                    readfile = False
                    if fn_match is not None:
                        for m in fn_match:
                            if m in fn:
                                readfile = True
                                break
                    else:
                        readfile = True
                        # Limit search to files under 1MB to avoid 30GB FastQ files etc.
                        try:
                            filesize = os.path.getsize(os.path.join(root,fn))
                        except (IOError, OSError, ValueError, UnicodeDecodeError):
                            log.debug("Couldn't read file when checking filesize: {}".format(fn))
                            readfile = False
                        else:
                            if filesize > 1000000:
                                readfile = False
                        
                        # Use mimetypes to exclude binary files where possible
                        (ftype, encoding) = mimetypes.guess_type(os.path.join(root, fn))
                        if encoding is not None:
                            readfile = False # eg. gzipped files
                        if ftype is not None and ftype.startswith('text') is False:
                            readfile = False # eg. images - 'image/jpeg'
                                
                    if readfile:
                        try:
                            with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                                # Search this file for our string of interest
                                returnfile = False
                                if contents_match is not None:
                                    for line in f:
                                        for m in contents_match:
                                            if m in line:
                                                returnfile = True
                                                break
                                    f.seek(0)
                                else:
                                    returnfile = True
                                
                                # Give back what was asked for. Yield instead of return
                                # so that this function can be used as an interator
                                # without loading all files at once.
                                if returnfile:
                                    if filehandles:
                                        yield {'s_name': s_name, 'f': f, 'root': root, 'fn': fn}
                                    else:
                                        yield {'s_name': s_name, 'f': f.read(), 'root': root, 'fn': fn}
                        except (IOError, OSError, ValueError, UnicodeDecodeError):
                            logger.debug("Couldn't read file when looking for output: {}".format(fn))
    
    
    def clean_s_name(self, s_name, root):
        """ Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary.
        :param s_name: The sample name to clean
        :param root: The directory path that this file is within
        :param prepend_dirs: boolean, whether to prepend dir name to s_name
        :param trimmed: boolean, remove common trimming suffixes from name?
        :return: The cleaned sample name, ready to be used
        """
        # Split then take first section to remove everything after these matches
        for ext in config.fn_clean_exts:
            s_name = s_name.split(ext ,1)[0]
        if config.prepend_dirs:
            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
        return s_name
    
    
    def general_stats_addcols(self, data, headers=None):
        """ Helper function to add to the General Statistics table.
        Adds to config.general_stats and does not return anything.
        :param data: A dict with the data. First key should be sample name,
                     then the data key, then the data.
        :param headers: Dict / OrderedDict with information for the headers, 
                        such as colour scales, min and max values etc.
                        See docs/writing_python.md for more information.
        :return: None
        """
        keys = data.keys()
        if headers is not None:
            keys = headers.keys()
        for k in keys:
            # Unique id to avoid overwriting by other modules
            if self.name is None:
                rid = '{}_{}'.format(''.join(random.sample(letters, 4)), k)
            else:
                safe_name = ''.join(c for c in self.name if c.isalnum()).lower()
                rid = '{}_{}'.format(safe_name, k)
            
            # Build the header cell
            try: scale = headers[k]['scale']
            except KeyError: scale = 'GnBu'
            
            try: dmax = 'data-chroma-max="{0}"'.format(headers[k]['max'])
            except KeyError: dmax = ''
            
            try: dmin = 'data-chroma-min="{0}"'.format(headers[k]['min'])
            except KeyError: dmin = ''
            
            try: title = headers[k]['title']
            except KeyError: title = k
            
            try: title = '<span data-toggle="tooltip" title="{}: {}">{}</span>'.format(self.name, headers[k]['description'], title)
            except KeyError: pass
            
            config.general_stats['headers'][rid] = '<th id="header_{}" class="chroma-col" data-chroma-scale="{}" {} {}>{}</th>'.format(rid, scale, dmax, dmin, title)
            
            # Add the data cells
            nrows = 0
            for (sname, samp) in data.items():
                if k in samp:
                    val = samp[k]
                    if 'modify' in headers[k] and callable(headers[k]['modify']):
                        val = headers[k]['modify'](val)
                        
                    formatstring = headers[k].get('')
                    try: formatstring = headers[k]['format']
                    except: formatstring = '{:.1f}'
                    try: val = formatstring.format(val)
                    except ValueError: val = formatstring.format(float(samp[k]))
                    except: val = samp[k]
                    
                    config.general_stats['rows'][sname][rid] = '<td class="{}">{}</td>'.format(rid, val)
                    nrows += 1
            
            # Remove header if we don't have any filled cells for it
            if nrows == 0:
                config.general_stats['headers'].pop(rid, None)
                logger.debug('Removing header {} from general stats table, as no data'.format(k))
        
        return None # it's good to be explicit, right?
        
        
    
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
                maxval = 0
                for k in sorted(d[s].keys()):
                    pairs.append([k, d[s][k]])
                    maxval = max(maxval, d[s][k])
                if maxval > 0 or config.get('hide_empty') is not True:
                    this_series = { 'name': s, 'data': pairs }
                    try:
                        this_series['color'] = config['colors'][s]
                    except: pass
                    thisplotdata.append(this_series)
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
            config['tt_label'] = 'Click to show original plot.<br>{}'.format(config.get('tt_label', '{point.x}: {point.y:.2f}'))
            if len(original_plots) > 1:
                next_prev_buttons = '<div class="clearfix"><div class="btn-group btn-group-sm"> \n\
                    <a href="#{prev}" class="btn btn-default original_plot_prev_btn" data-target="#{id}">&laquo; Previous</a> \n\
                    <a href="#{next}" class="btn btn-default original_plot_nxt_btn" data-target="#{id}">Next &raquo;</a> \n\
                </div></div>'.format(id=config['id'], prev=original_plots[-1]['s_name'], next=original_plots[1]['s_name'])
            else:
                next_prev_buttons = ''
            html += '<p class="text-muted instr">Click to show original FastQC plot.</p>\n\
                    <div id="{id}_wrapper" class="hc-plot-wrapper"> \n\
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
            html += '<div id="{id}" class="hc-plot hc-line-plot"></div> \n'.format(id=config['id'])
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
                     Can supply a list of dicts and will have buttons to switch
        :param cats: optional list, dict or OrderedDict with plot categories
        :param config: optional dict with config key:value pairs
        :return: HTML and JS, ready to be inserted into the page
        """
        
        # Given one dataset - turn it into a list
        if type(data) is not list:
            data = [data]
        
        # Check we have a list of cats
        if type(cats) is not list or type(cats[0]) is str:
            cats = [cats]
        
        # Check that we have cats at all - find them from the data
        for idx, cat in enumerate(cats):
            if cats[idx] is None:
                cats[idx] = list(set(k for s in data[idx].keys() for k in data[idx][s].keys() ))
        
        # Given a list of cats - turn it into a dict
        for idx, cat in enumerate(cats):
            if type(cat) is list:
                newcats = OrderedDict()
                for c in cat:
                    newcats[c] = {'name': c}
                cats[idx] = newcats
        
        # Parse the data into a HighCharts friendly format
        plotsamples = list()
        plotdata = list()
        for idx, d in enumerate(data):
            hc_samples = sorted(list(d.keys()))
            hc_data = list()
            for c in cats[idx].keys():
                thisdata = list()
                for s in hc_samples:
                    thisdata.append(d[s][c])
                if max(thisdata) > 0:
                    thisdict = { 'name': cats[idx][c]['name'], 'data': thisdata }
                    if 'color' in cats[idx][c]:
                        thisdict['color'] = cats[idx][c]['color']
                    hc_data.append(thisdict)
            plotsamples.append(hc_samples)
            plotdata.append(hc_data)
        
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
    		</div> '.format(id=config['id'], c_a=c_active, p_a=p_active, c_l=c_label, p_l=p_label)
            if len(plotdata) > 1:
                html += ' &nbsp; &nbsp; '
        
        # Buttons to cycle through different datasets
        if len(plotdata) > 1:
            html += '<div class="btn-group switch_group">\n'
            for k, p in enumerate(plotdata):
                active = 'active' if k == 0 else ''
                try: name = config['data_labels'][k]
                except: name = k+1
                try: ylab = 'data-ylab="{}"'.format(config['data_labels'][k]['ylab'])
                except: ylab = 'data-ylab="{}"'.format(name) if name != k+1 else ''
                html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} data-newdata="{id}_datasets[{k}]" data-target="#{id}">{n}</button>\n'.format(a=active, id=config['id'], n=name, y=ylab, k=k)
            html += '</div>\n\n'
        
        # Plot and javascript function
        html += '<div id="{id}" class="hc-plot hc-bar-plot"></div> \n\
        <script type="text/javascript"> \n\
            var {id}_samples = {s}; \n\
            var {id}_datasets = {d}; \n\
            $(function () {{ plot_stacked_bar_graph("#{id}", {id}_samples[0], {id}_datasets[0], {c}); }}); \
        </script>'.format(id=config['id'], s=json.dumps(plotsamples), d=json.dumps(plotdata), c=json.dumps(config));
        
        return html
        
    
    def write_csv_file(self, data, fn):
        with io.open (os.path.join(config.output_dir, 'report_data', fn), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( data ), file=f)
            
            
    def dict_to_csv (self, d, delim="\t"):
        """ Converts a dict to a CSV string
        :param d: 2D dictionary, first keys sample names and second key
                  column headers
        :param delim: optional delimiter character. Default: \t
        :return: Flattened string, suitable to write to a CSV file.
        """

        h = None # We make a list of keys to ensure consistent order
        l = list()
        for sn in sorted(d.keys()):
            if h is None:
                h = list(d[sn].keys())
                l.append(delim.join([''] + h))
            thesefields = [sn] + [ str(d[sn].get(k, '')) for k in h ]
            l.append( delim.join( thesefields ) )
        return ('\n'.join(l)).encode('utf-8', 'ignore').decode('utf-8')
