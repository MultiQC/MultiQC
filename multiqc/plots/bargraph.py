#!/usr/bin/env python

""" MultiQC functions to plot a bargraph """

from __future__ import print_function
import base64
from collections import OrderedDict
import io
import json
import logging
import math
import os
import random
import sys

from multiqc.utils import config
logger = logging.getLogger(__name__)

try:
    # Import matplot lib but avoid default X environment
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception as e:
    # MatPlotLib can break in a variety of ways. Fake an error message and continue without it if so.
    # The lack of the library will be handled when plots are attempted
    print("##### ERROR! MatPlotLib library could not be loaded!    #####", file=sys.stderr)
    print("##### Flat plots will instead be plotted as interactive #####", file=sys.stderr)
    logger.exception(e)

letters = 'abcdefghijklmnopqrstuvwxyz'

# Load the template so that we can access it's configuration
template_mod = config.avail_templates[config.template].load()

def plot (data, cats=None, pconfig={}):
    """ Plot a horizontal bar graph. Expects a 2D dict of sample
    data. Also can take info about categories. There are quite a
    few variants of how to use this function, see the docs for details.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
                 Can supply a list of dicts and will have buttons to switch
    :param cats: optional list, dict or OrderedDict with plot categories
    :param pconfig: optional dict with config key:value pairs
    :return: HTML and JS, ready to be inserted into the page
    """
    
    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]
    
    # Check we have a list of cats
    try:
        cats[0].keys()
    except (KeyError, AttributeError):
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
    
    # Parse the data into a chart friendly format
    plotsamples = list()
    plotdata = list()
    for idx, d in enumerate(data):
        hc_samples = sorted(list(d.keys()))
        hc_data = list()
        sample_dcount = dict()
        for c in cats[idx].keys():
            thisdata = list()
            for s in hc_samples:
                try:
                    thisdata.append(d[s][c])
                    try:
                        sample_dcount[s] += 1
                    except KeyError:
                        sample_dcount[s] = 1
                except KeyError:
                    pass
            if len(thisdata) > 0 and max(thisdata) > 0:
                thisdict = { 'name': cats[idx][c]['name'], 'data': thisdata }
                if 'color' in cats[idx][c]:
                    thisdict['color'] = cats[idx][c]['color']
                hc_data.append(thisdict)
        if len(hc_data) > 0:
            plotsamples.append([s for s in hc_samples if sample_dcount.get(s, 0) > 0])
            plotdata.append(hc_data)
    
    if len(plotdata) == 0:
        logger.warning('Tried to make bar plot, but had no data')
        return '<p class="text-danger">Error - was not able to plot data.</p>'
    
    # Make a plot - custom, interactive or flat
    try:
        return template_mod.bargraph(plotdata, plotsamples, pconfig)
    except (AttributeError, TypeError):
        if config.plots_force_flat or (not config.plots_force_interactive and len(plotsamples[0]) > config.plots_flat_numseries):
            try:
                return matplotlib_bargraph(plotdata, plotsamples, pconfig)
            except:
                logger.error("############### Error making MatPlotLib figure! Falling back to HighCharts.")
                return highcharts_bargraph(plotdata, plotsamples, pconfig)
        else:
            return highcharts_bargraph(plotdata, plotsamples, pconfig)



def highcharts_bargraph (plotdata, plotsamples=None, pconfig={}):
    """
    Build the HTML needed for a HighCharts bar graph. Should be
    called by plot_bargraph, which properly formats input data.
    """
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_hcplot_'+''.join(random.sample(letters, 10))
    html = '<div class="mqc_hcplot_plotgroup">'
    
    # Counts / Percentages / Log Switches
    if pconfig.get('cpswitch') is not False or pconfig.get('logswitch') is True:
        if pconfig.get('cpswitch_c_active', True) is True:
            c_active = 'active'
            p_active = ''
            l_active = ''
        elif pconfig.get('logswitch_active') is True:
            c_active = ''
            p_active = ''
            l_active = 'active'
        else:
            c_active = ''
            p_active = 'active'
            l_active = ''
            pconfig['stacking'] = 'percent'
        c_label = pconfig.get('cpswitch_counts_label', 'Counts')
        p_label = pconfig.get('cpswitch_percent_label', 'Percentages')
        l_label = pconfig.get('logswitch_label', 'Log10')
        html += '<div class="btn-group hc_switch_group"> \n'
        html += '<button class="btn btn-default btn-sm {c_a}" data-action="set_numbers" data-target="{id}" data-ylab="{c_l}">{c_l}</button> \n'.format(id=pconfig['id'], c_a=c_active, c_l=c_label)
        if pconfig.get('cpswitch', True) is True:
            html += '<button class="btn btn-default btn-sm {p_a}" data-action="set_percent" data-target="{id}" data-ylab="{p_l}">{p_l}</button> \n'.format(id=pconfig['id'], p_a=p_active, p_l=p_label)
        if pconfig.get('logswitch') is True:
            html += '<button class="btn btn-default btn-sm {l_a}" data-action="set_log" data-target="{id}" data-ylab="{l_l}">{l_l}</button> \n'.format(id=pconfig['id'], l_a=l_active, l_l=l_label)
            pconfig['reversedStacks'] = True
        html += '</div> '
        if len(plotdata) > 1:
            html += ' &nbsp; &nbsp; '
    
    # Buttons to cycle through different datasets
    if len(plotdata) > 1:
        html += '<div class="btn-group hc_switch_group">\n'
        for k, p in enumerate(plotdata):
            active = 'active' if k == 0 else ''
            try: name = pconfig['data_labels'][k]
            except: name = k+1
            try: ylab = 'data-ylab="{}"'.format(pconfig['data_labels'][k]['ylab'])
            except: ylab = 'data-ylab="{}"'.format(name) if name != k+1 else ''
            try: ymax = 'data-ymax="{}"'.format(pconfig['data_labels'][k]['ymax'])
            except: ymax = ''
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(a=active, id=pconfig['id'], n=name, y=ylab, k=k)
        html += '</div>\n\n'
    
    # Plot and javascript function
    html += '<div class="hc-plot-wrapper"><div id="{id}" class="hc-plot not_rendered hc-bar-plot"><small>loading..</small></div></div> \n\
    </div> \n\
    <script type="text/javascript"> \n\
        mqc_plots["{id}"] = {{ \n\
            "plot_type": "bar_graph", \n\
            "samples": {s}, \n\
            "datasets": {d}, \n\
            "config": {c} \n\
        }} \n\
    </script>'.format(id=pconfig['id'], s=json.dumps(plotsamples), d=json.dumps(plotdata), c=json.dumps(pconfig));
    
    return html


def matplotlib_bargraph (plotdata, plotsamples, pconfig={}):
    """
    Plot a bargraph with Matplot lib and return a HTML string. Either embeds a base64
    encoded image within HTML or writes the plot and links to it. Should be called by
    plot_bargraph, which properly formats the input data.
    """
    
    # Plot group ID
    if pconfig.get('id') is None:
        pconfig['id'] = 'mqc_mplplot_'+''.join(random.sample(letters, 10))
    # Individual plot IDs
    pids = []
    for k in range(len(plotdata)):
        try:
            name = pconfig['data_labels'][k]
        except:
            name = k+1
        pid = 'mqc_{}_{}'.format(pconfig['id'], name)
        pid = "".join([c for c in pid if c.isalpha() or c.isdigit() or c == '_' or c == '-'])
        pids.append(pid)
    
    html = '<div class="mqc_mplplot_plotgroup" id="{}">'.format(pconfig['id'])
    
    # Same defaults as HighCharts for consistency
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9', 
                      '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']
    
    # Counts / Percentages Switch
    if pconfig.get('cpswitch') is not False:
        if pconfig.get('cpswitch_c_active', True) is True:
            c_active = 'active'
            p_active = ''
        else:
            c_active = ''
            p_active = 'active'
            pconfig['stacking'] = 'percent'
        c_label = pconfig.get('cpswitch_counts_label', 'Counts')
        p_label = pconfig.get('cpswitch_percent_label', 'Percentages')
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_setcountspcnt"> \n\
            <button class="btn btn-default btn-sm {c_a} counts">{c_l}</button> \n\
            <button class="btn btn-default btn-sm {p_a} pcnt">{p_l}</button> \n\
        </div> '.format(c_a=c_active, p_a=p_active, c_l=c_label, p_l=p_label)
        if len(plotdata) > 1:
            html += ' &nbsp; &nbsp; '
    
    # Buttons to cycle through different datasets
    if len(plotdata) > 1:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(plotdata):
            pid = pids[k]
            active = 'active' if k == 0 else ''
            try:
                name = pconfig['data_labels'][k]
            except:
                name = k+1
            html += '<button class="btn btn-default btn-sm {a}" data-target="#{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name)
        html += '</div>\n\n'
    
    # Go through datasets creating plots
    for pidx, pdata in enumerate(plotdata):
        
        # Plot percentage as well as counts
        plot_pcts = [False]
        if pconfig.get('cpswitch') is not False:
            plot_pcts = [False, True]
        
        for plot_pct in plot_pcts:
            
            # Plot ID
            pid = pids[pidx]
            hide_plot = False
            if plot_pct is True:
                pid = '{}_pc'.format(pid)
                if pconfig.get('cpswitch_c_active', True) is True:
                    hide_plot = True
            else:
                if pconfig.get('cpswitch_c_active', True) is not True:
                    hide_plot = True
            
            # Set up figure
            plt_height = min(30, max(6, len(plotsamples[pidx]) / 2.3))
            fig = plt.figure(figsize=(14, plt_height), frameon=False)
            axes = fig.add_subplot(111)
            y_ind = range(len(plotsamples[pidx]))
            bar_width = 0.8
            
            # Count totals for each sample
            if plot_pct is True:
                s_totals = [0 for _ in pdata[0]['data']]
                for series_idx, d in enumerate(pdata):
                    for sample_idx, v in enumerate(d['data']):
                        s_totals[sample_idx] += v
            
            # Plot bars
            dlabels = []
            for idx, d in enumerate(pdata):
                
                # Plot percentages
                values = d['data']
                if len(values) < len(y_ind):
                    values.extend([0] * (len(y_ind) - len(values)))
                if plot_pct is True:
                    for (key,var) in enumerate(values):
                        s_total = s_totals[key]
                        if s_total == 0:
                            values[key] = 0
                        else:
                            values[key] = (float(var+0.0)/float(s_total))*100
                
                # Get offset for stacked bars
                if idx == 0:
                    prevdata = [0] * len(plotsamples[pidx])
                else:
                    for i, p in enumerate(prevdata):
                        prevdata[i] += pdata[idx-1]['data'][i]
                # Default colour index
                cidx = idx
                while cidx >= len(default_colors):
                    cidx -= len(default_colors)
                # Save the name of this series
                dlabels.append(d['name'])
                # Add the series of bars to the plot
                axes.barh(
                    y_ind, values, bar_width, left=prevdata,
                    color=d.get('color', default_colors[cidx]), align='center', linewidth=1, edgecolor='w'
                )
            
            # Tidy up axes
            axes.tick_params(labelsize=8, direction='out', left=False, right=False, top=False, bottom=False)
            axes.set_xlabel(pconfig.get('ylab', '')) # I know, I should fix the fact that the config is switched
            axes.set_ylabel(pconfig.get('xlab', ''))
            axes.set_yticks(y_ind) # Specify where to put the labels
            axes.set_yticklabels(plotsamples[pidx]) # Set y axis sample name labels
            axes.set_ylim((-0.5, len(y_ind)-0.5)) # Reduce padding around plot area
            if plot_pct is True:
                axes.set_xlim((0, 100))
                # Add percent symbols
                vals = axes.get_xticks()
                axes.set_xticklabels(['{:.0f}%'.format(x) for x in vals])
            else:
                default_xlimits = axes.get_xlim()
                axes.set_xlim((pconfig.get('ymin', default_xlimits[0]),pconfig.get('ymax', default_xlimits[1])))
            if 'title' in pconfig:
                top_gap = 1 + (0.5 / plt_height)
                plt.text(0.5, top_gap, pconfig['title'], horizontalalignment='center', fontsize=16, transform=axes.transAxes)
            axes.grid(True, zorder=0, which='both', axis='x', linestyle='-', color='#dedede', linewidth=1)
            axes.set_axisbelow(True)
            axes.spines['right'].set_visible(False)
            axes.spines['top'].set_visible(False)
            axes.spines['bottom'].set_visible(False)
            axes.spines['left'].set_visible(False)
            plt.gca().invert_yaxis() # y axis is reverse sorted otherwise
            
            # Hide some labels if we have a lot of samples
            show_nth = max(1, math.ceil(len(pdata[0]['data'])/150))
            for idx, label in enumerate(axes.get_yticklabels()):
                if idx % show_nth != 0:
                    label.set_visible(False)
            
            # Legend
            bottom_gap = -1 * (1 - ((plt_height - 1.5) / plt_height))
            lgd = axes.legend(dlabels, loc='lower center', bbox_to_anchor=(0, bottom_gap, 1, .102), ncol=5, mode='expand', fontsize=8, frameon=False)
            
            # Should this plot be hidden on report load?
            hidediv = ''
            if pidx > 0 or hide_plot:
                hidediv = ' style="display:none;"'
            
            # Save the plot to the data directory
            plot_dir = os.path.join(config.data_dir, 'multiqc_plots')
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
            plot_fn = os.path.join(plot_dir, '{}.png'.format(pid))
            fig.savefig(plot_fn, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
            
            # Output the figure to a base64 encoded string
            if getattr(template_mod, 'base64_plots', True) is True:
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format='png', bbox_inches='tight')
                b64_img = base64.b64encode(img_buffer.getvalue()).decode('utf8')
                img_buffer.close()
                html += '<div class="mqc_mplplot" id="{}"{}><img src="data:image/png;base64,{}" /></div>'.format(pid, hidediv, b64_img)
            
            # Link to the saved image
            else:
                plot_relpath = os.path.join(config.data_dir_name, 'multiqc_plots', '{}.png'.format(pid))
                html += '<div class="mqc_mplplot" id="{}"{}><img src="{}" /></div>'.format(pid, hidediv, plot_relpath)
            
            plt.close(fig)
            
    
    # Close wrapping div
    html += '</div>'
    
    return html

