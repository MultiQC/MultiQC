#!/usr/bin/env python

""" MultiQC functions to plot a linegraph """

from __future__ import print_function
from collections import OrderedDict
import base64
import io
import json
import logging
import os
import random
import sys

from multiqc.utils import config, report, util_functions
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

def plot (data, pconfig={}):
    """ Plot a line graph with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """

    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]

    # Smooth dataset if requested in config
    if pconfig.get('smooth_points', None) is not None:
        sumcounts = pconfig.get('smooth_points_sumcounts', True)
        for i, d in enumerate(data):
            sumc = sumcounts
            if type(sumcounts) is list:
                sumc = sumcounts[i]
            data[i] = smooth_line_data(d, pconfig['smooth_points'], sumc)

    # Generate the data dict structure expected by HighCharts series
    plotdata = list()
    for d in data:
        thisplotdata = list()
        for s in sorted(d.keys()):
            pairs = list()
            maxval = 0
            if 'categories' in pconfig:
                pconfig['categories'] = list()
                for k in d[s].keys():
                    pconfig['categories'].append(k)
                    pairs.append(d[s][k])
                    maxval = max(maxval, d[s][k])
            else:
                for k in sorted(d[s].keys()):
                    pairs.append([k, d[s][k]])
                    try:
                        maxval = max(maxval, d[s][k])
                    except TypeError:
                        pass
            if maxval > 0 or pconfig.get('hide_empty') is not True:
                this_series = { 'name': s, 'data': pairs }
                try:
                    this_series['color'] = pconfig['colors'][s]
                except: pass
                thisplotdata.append(this_series)
        plotdata.append(thisplotdata)

    # Add on annotation data series
    try:
        if 'extra_series' in pconfig:
            extra_series = pconfig['extra_series']
            if type(pconfig['extra_series']) == dict:
                extra_series = [[ pconfig['extra_series'] ]]
            elif type(pconfig['extra_series']) == list and type(pconfig['extra_series'][0]) == dict:
                extra_series = [ pconfig['extra_series'] ]
            for i, es in enumerate(extra_series):
                for s in es:
                    plotdata[i].append(s)
    except KeyError:
        pass

    # Make a plot - template custom, or interactive or flat
    try:
        return template_mod.linegraph(plotdata, pconfig)
    except (AttributeError, TypeError):
        if config.plots_force_flat or (not config.plots_force_interactive and len(plotdata[0]) > config.plots_flat_numseries):
            try:
                return matplotlib_linegraph(plotdata, pconfig)
            except:
                logger.error("############### Error making MatPlotLib figure! Falling back to HighCharts.")
                return highcharts_linegraph(plotdata, pconfig)
        else:
            # Use MatPlotLib to generate static plots if requested
            if config.export_plots:
                matplotlib_linegraph(plotdata, pconfig)
            # Return HTML for HighCharts dynamic plot
            return highcharts_linegraph(plotdata, pconfig)



def highcharts_linegraph (plotdata, pconfig={}):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by linegraph.plot(), which properly formats input data.
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
            try:
                name = pconfig['data_labels'][k]['name']
            except:
                name = k+1
            try:
                ylab = 'data-ylab="{}"'.format(pconfig['data_labels'][k]['ylab'])
            except:
                ylab = 'data-ylab="{}"'.format(name) if name != k+1 else ''
            try:
                ymax = 'data-ymax="{}"'.format(pconfig['data_labels'][k]['ymax'])
            except:
                ymax = ''
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(a=active, id=pconfig['id'], n=name, y=ylab, ym=ymax, k=k)
        html += '</div>\n\n'

    # The plot div
    html += '<div class="hc-plot-wrapper"><div id="{id}" class="hc-plot not_rendered hc-line-plot"><small>loading..</small></div></div></div> \n'.format(id=pconfig['id'])

    # Javascript with data dump
    html += '<script type="text/javascript"> \n\
        mqc_plots["{id}"] = {{ \n\
            "plot_type": "xy_line", \n\
            "datasets": {d}, \n\
            "config": {c} \n\
        }} \n\
    </script>'.format(id=pconfig['id'], d=json.dumps(plotdata), c=json.dumps(pconfig));

    report.num_hc_plots += 1
    return html


def matplotlib_linegraph (plotdata, pconfig={}):
    """
    Plot a line graph with Matplot lib and return a HTML string. Either embeds a base64
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
            name = pconfig['data_labels'][k]['name']
        except:
            name = k+1
        pid = 'mqc_{}_{}'.format(pconfig['id'], name)
        pid = "".join([c for c in pid if c.isalpha() or c.isdigit() or c == '_' or c == '-'])
        pids.append(pid)

    html = '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> ' + \
          'Flat image plot. Toolbox functions such as highlighting / hiding samples will not work ' + \
          '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    html += '<div class="mqc_mplplot_plotgroup" id="{}">'.format(pconfig['id'])

    # Same defaults as HighCharts for consistency
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
                      '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']

    # Buttons to cycle through different datasets
    if len(plotdata) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(plotdata):
            pid = pids[k]
            active = 'active' if k == 0 else ''
            try:
                name = pconfig['data_labels'][k]['name']
            except:
                name = k+1
            html += '<button class="btn btn-default btn-sm {a}" data-target="#{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name)
        html += '</div>\n\n'

    # Go through datasets creating plots
    for pidx, pdata in enumerate(plotdata):

        # Plot ID
        pid = pids[pidx]

        # Save plot data to file
        fdata = OrderedDict()
        lastcats = None
        sharedcats = True
        for d in pdata:
            fdata[d['name']] = OrderedDict()
            for i, x in enumerate(d['data']):
                if type(x) is list:
                    fdata[d['name']][str(x[0])] = x[1]
                    # Check to see if all categories are the same
                    if lastcats is None:
                        lastcats = [x[0] for x in d['data']]
                    elif lastcats != [x[0] for x in d['data']]:
                        sharedcats = False
                else:
                    try:
                        fdata[d['name']][pconfig['categories'][i]] = x
                    except (KeyError, IndexError):
                        fdata[d['name']][str(i)] = x
        # Custom tsv output if the x axis varies
        if not sharedcats and config.data_format == 'tsv':
            fout = ''
            for d in pdata:
                fout += "\t"+"\t".join([str(x[0]) for x in d['data']])
                fout += "\n{}\t".format(d['name'])
                fout += "\t".join([str(x[1]) for x in d['data']])
                fout += "\n"
            with io.open (os.path.join(config.data_dir, '{}.txt'.format(pid)), 'w', encoding='utf-8') as f:
                print( fout.encode('utf-8', 'ignore').decode('utf-8'), file=f )
        else:
            util_functions.write_data_file(fdata, pid)

        # Set up figure
        fig = plt.figure(figsize=(14, 6), frameon=False)
        axes = fig.add_subplot(111)

        # Go through data series
        for idx, d in enumerate(pdata):

            # Default colour index
            cidx = idx
            while cidx >= len(default_colors):
                cidx -= len(default_colors)

            # Line style
            linestyle = 'solid'
            if d.get('dashStyle', None) == 'Dash':
                linestyle = 'dashed'

            # Reformat data (again)
            try:
                axes.plot([x[0] for x in d['data']], [x[1] for x in d['data']], label=d['name'], color=d.get('color', default_colors[cidx]), linestyle=linestyle, linewidth=1, marker=None)
            except TypeError:
                # Categorical data on x axis
                axes.plot(d['data'], label=d['name'], color=d.get('color', default_colors[cidx]), linewidth=1, marker=None)

        # Tidy up axes
        axes.tick_params(labelsize=8, direction='out', left=False, right=False, top=False, bottom=False)
        axes.set_xlabel(pconfig.get('xlab', ''))
        axes.set_ylabel(pconfig.get('ylab', ''))

        # Dataset specific y label
        try:
            axes.set_ylabel(pconfig['data_labels'][pidx]['ylab'])
        except:
            pass

        # Axis limits
        default_ylimits = axes.get_ylim()
        ymin = default_ylimits[0]
        if 'ymin' in pconfig:
            ymin = pconfig['ymin']
        elif 'yCeiling' in pconfig:
            ymin = min(pconfig['yCeiling'], default_ylimits[0])
        ymax = default_ylimits[1]
        if 'ymax' in pconfig:
            ymax = pconfig['ymax']
        elif 'yFloor' in pconfig:
            ymax = max(pconfig['yCeiling'], default_ylimits[1])
        if (ymax - ymin) < pconfig.get('yMinRange', 0):
            ymax = ymin + pconfig['yMinRange']
        axes.set_ylim((ymin, ymax))

        # Dataset specific ymax
        try:
            axes.set_ylim((ymin, pconfig['data_labels'][pidx]['ymax']))
        except:
            pass

        default_xlimits = axes.get_xlim()
        xmin = default_xlimits[0]
        if 'xmin' in pconfig:
            xmin = pconfig['xmin']
        elif 'xCeiling' in pconfig:
            xmin = min(pconfig['xCeiling'], default_xlimits[0])
        xmax = default_xlimits[1]
        if 'xmax' in pconfig:
            xmax = pconfig['xmax']
        elif 'xFloor' in pconfig:
            xmax = max(pconfig['xCeiling'], default_xlimits[1])
        if (xmax - xmin) < pconfig.get('xMinRange', 0):
            xmax = xmin + pconfig['xMinRange']
        axes.set_xlim((xmin, xmax))

        # Plot title
        if 'title' in pconfig:
            plt.text(0.5, 1.05, pconfig['title'], horizontalalignment='center', fontsize=16, transform=axes.transAxes)
        axes.grid(True, zorder=10, which='both', axis='y', linestyle='-', color='#dedede', linewidth=1)

        # X axis categories, if specified
        if 'categories' in pconfig:
            axes.set_xticks([i for i,v in enumerate(pconfig['categories'])])
            axes.set_xticklabels(pconfig['categories'])

        # Axis lines
        xlim = axes.get_xlim()
        axes.plot([xlim[0], xlim[1]], [0, 0], linestyle='-', color='#dedede', linewidth=2)
        axes.set_axisbelow(True)
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['bottom'].set_visible(False)
        axes.spines['left'].set_visible(False)

        # Background colours, if specified
        if 'yPlotBands' in pconfig:
            xlim = axes.get_xlim()
            for pb in pconfig['yPlotBands']:
                axes.barh(pb['from'], xlim[1], height = pb['to']-pb['from'], left=xlim[0], color=pb['color'], linewidth=0, zorder=0)
        if 'xPlotBands' in pconfig:
            ylim = axes.get_ylim()
            for pb in pconfig['xPlotBands']:
                axes.bar(pb['from'], ylim[1], width = pb['to']-pb['from'], bottom=ylim[0], color=pb['color'], linewidth=0, zorder=0)

        # Tight layout - makes sure that legend fits in and stuff
        if len(pdata) <= 15:
            lgd = axes.legend(loc='lower center', bbox_to_anchor=(0, -0.22, 1, .102), ncol=5, mode='expand', fontsize=8, frameon=False)
            plt.tight_layout(rect=[0,0.08,1,0.92])
        else:
            plt.tight_layout(rect=[0,0,1,0.92])

        # Should this plot be hidden on report load?
        hidediv = ''
        if pidx > 0:
            hidediv = ' style="display:none;"'

        # Save the plot to the data directory if export is requests
        if config.export_plots:
            for fformat in config.export_plot_formats:
                # Make the directory if it doesn't already exist
                plot_dir = os.path.join(config.plots_dir, fformat)
                if not os.path.exists(plot_dir):
                    os.makedirs(plot_dir)
                # Save the plot
                plot_fn = os.path.join(plot_dir, '{}.{}'.format(pid, fformat))
                fig.savefig(plot_fn, format=fformat, bbox_inches='tight')

        # Output the figure to a base64 encoded string
        if getattr(template_mod, 'base64_plots', True) is True:
            img_buffer = io.BytesIO()
            fig.savefig(img_buffer, format='png', bbox_inches='tight')
            b64_img = base64.b64encode(img_buffer.getvalue()).decode('utf8')
            img_buffer.close()
            html += '<div class="mqc_mplplot" id="{}"{}><img src="data:image/png;base64,{}" /></div>'.format(pid, hidediv, b64_img)

        # Save to a file and link <img>
        else:
            plot_relpath = os.path.join(config.data_dir_name, 'multiqc_plots', '{}.png'.format(pid))
            html += '<div class="mqc_mplplot" id="{}"{}><img src="{}" /></div>'.format(pid, hidediv, plot_relpath)

        plt.close(fig)


    # Close wrapping div
    html += '</div>'

    report.num_mpl_plots += 1
    return html



def smooth_line_data(data, numpoints, sumcounts=True):
    """
    Function to take an x-y dataset and use binning to
    smooth to a maximum number of datapoints.
    """
    smoothed = {}
    for s_name, d in data.items():
        smoothed[s_name] = OrderedDict();
        p = 0
        binsize = len(d) / numpoints
        if binsize < 1:
            binsize = 1
        binvals = []
        for x, y in d.items():
            if p < binsize:
                binvals.append(y)
                p += 1
            else:
                if sumcounts is True:
                    v = sum(binvals)
                else:
                    v = sum(binvals) / binsize
                smoothed[s_name][x] = v
                p = 0
                binvals = []
    return smoothed

