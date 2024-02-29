#!/usr/bin/env python

""" MultiQC functions to plot a bargraph """


import base64
import io
import logging

import math
import os
import random
import sys

from multiqc.utils import config, report, util_functions

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(plotdata, plotsamples, pconfig=None):
    """Plot a horizontal bar graph. Expects a 2D dict of sample
    data. Also, can take info about categories. There are quite a
    few variants of how to use this function, see the docs for details.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
                 Can supply a list of dicts and will have buttons to switch
    :param cats: optional list or dict with plot categories
    :param pconfig: optional dict with config key:value pairs
    :return: HTML and JS, ready to be inserted into the page
    """
    if config.plots_force_flat or (
        not config.plots_force_interactive and len(plotsamples[0]) > config.plots_flat_numseries
    ):
        try:
            report.num_mpl_plots += 1
            return matplotlib_bargraph(plotdata, plotsamples, pconfig)
        except Exception as e:
            if config.strict:
                raise
            logger.error("############### Error making MatPlotLib figure! Falling back to HighCharts.")
            logger.debug(e, exc_info=True)
            return highcharts_bargraph(plotdata, plotsamples, pconfig)
    else:
        # Use MatPlotLib to generate static plots if requested
        if config.export_plots:
            try:
                matplotlib_bargraph(plotdata, plotsamples, pconfig)
            except Exception as e:
                logger.error("############### Error making MatPlotLib figure! Plot not exported.")
                logger.debug(e, exc_info=True)
        # Return HTML for HighCharts dynamic plot
        return highcharts_bargraph(plotdata, plotsamples, pconfig)


def highcharts_bargraph(plotdata, plotsamples=None, pconfig=None):
    """
    Build the HTML needed for a HighCharts bar graph. Should be
    called by plot_bargraph, which properly formats input data.
    """
    if pconfig is None:
        pconfig = {}
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_hcplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    html = '<div class="mqc_hcplot_plotgroup">'

    # Counts / Percentages / Log Switches
    if pconfig.get("cpswitch") is not False or pconfig.get("logswitch") is True:
        if pconfig.get("logswitch_active") is True:
            c_active = ""
            p_active = ""
            l_active = "active"
        elif pconfig.get("cpswitch_c_active", True) is True:
            c_active = "active"
            p_active = ""
            l_active = ""
        else:
            c_active = ""
            p_active = "active"
            l_active = ""
            pconfig["stacking"] = "percent"
        c_label = pconfig.get("cpswitch_counts_label", "Counts")
        p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        l_label = pconfig.get("logswitch_label", "Log10")
        html += '<div class="btn-group hc_switch_group"> \n'
        html += '<button class="btn btn-default btn-sm {c_a}" data-action="set_numbers" data-target="{id}" data-ylab="{c_l}">{c_l}</button> \n'.format(
            id=pconfig["id"], c_a=c_active, c_l=c_label
        )
        if pconfig.get("cpswitch", True) is True:
            html += '<button class="btn btn-default btn-sm {p_a}" data-action="set_percent" data-target="{id}" data-ylab="{p_l}">{p_l}</button> \n'.format(
                id=pconfig["id"], p_a=p_active, p_l=p_label
            )
        if pconfig.get("logswitch") is True:
            html += '<button class="btn btn-default btn-sm {l_a}" data-action="set_log" data-target="{id}" data-ylab="{l_l}">{l_l}</button> \n'.format(
                id=pconfig["id"], l_a=l_active, l_l=l_label
            )
            pconfig["reversedStacks"] = True
        html += "</div> "
        if len(plotdata) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(plotdata) > 1:
        html += '<div class="btn-group hc_switch_group">\n'
        for k, p in enumerate(plotdata):
            active = "active" if k == 0 else ""
            try:
                name = pconfig["data_labels"][k]["name"]
            except Exception:
                try:
                    name = pconfig["data_labels"][k]
                except Exception:
                    name = k + 1
            try:
                ylab = f"data-ylab=\"{pconfig['data_labels'][k]['ylab']}\""
            except Exception:
                ylab = f'data-ylab="{name}"' if name != k + 1 else ""
            try:
                ymax = f"data-ymax=\"{pconfig['data_labels'][k]['ymax']}\""
            except Exception:
                ymax = ""
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(
                a=active, id=pconfig["id"], n=name, y=ylab, ym=ymax, k=k
            )
        html += "</div>\n\n"

    # Plot HTML
    html += """<div class="hc-plot-wrapper"{height}>
        <div id="{id}" class="hc-plot not_rendered hc-bar-plot"><small>loading..</small></div>
    </div></div>""".format(
        id=pconfig["id"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[pconfig["id"]] = {
        "plot_type": "bar_graph",
        "samples": plotsamples,
        "datasets": plotdata,
        "config": pconfig,
    }

    return html


def matplotlib_bargraph(plotdata, plotsamples, pconfig=None):
    """
    Plot a bargraph with Matplot lib and return a HTML string. Either embeds a base64
    encoded image within HTML or writes the plot and links to it. Should be called by
    plot_bargraph, which properly formats the input data.
    """
    try:
        # Import matplot lib but avoid default X environment
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        logger.debug(f"Using matplotlib version {matplotlib.__version__}")
    except Exception as e:
        # MatPlotLib can break in a variety of ways. Fake an error message and continue without it if so.
        # The lack of the library will be handled when plots are attempted
        print("##### ERROR! MatPlotLib library could not be loaded!    #####", file=sys.stderr)
        print("##### Flat plots will instead be plotted as interactive #####", file=sys.stderr)
        print(e)

    if pconfig is None:
        pconfig = {}

    # Plot group ID
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_mplplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    # Individual plot IDs
    pids = []
    for k in range(len(plotdata)):
        try:
            name = pconfig["data_labels"][k]
        except Exception:
            name = k + 1
        pid = f"mqc_{pconfig['id']}_{name}"
        pid = report.save_htmlid(pid, skiplint=True)
        pids.append(pid)

    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )
    html += f"<div class=\"mqc_mplplot_plotgroup\" id=\"{pconfig['id']}\">"

    # Counts / Percentages Switch
    if pconfig.get("cpswitch") is not False and not config.simple_output:
        if pconfig.get("cpswitch_c_active", True) is True:
            c_active = "active"
            p_active = ""
        else:
            c_active = ""
            p_active = "active"
            pconfig["stacking"] = "percent"
        c_label = pconfig.get("cpswitch_counts_label", "Counts")
        p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_setcountspcnt"> \n\
            <button class="btn btn-default btn-sm {c_a} counts">{c_l}</button> \n\
            <button class="btn btn-default btn-sm {p_a} pcnt">{p_l}</button> \n\
        </div> '.format(c_a=c_active, p_a=p_active, c_l=c_label, p_l=p_label)
        if len(plotdata) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(plotdata) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(plotdata):
            pid = pids[k]
            active = "active" if k == 0 else ""
            try:
                name = pconfig["data_labels"][k]
            except Exception:
                name = k + 1
            html += '<button class="btn btn-default btn-sm {a}" data-target="#{pid}">{n}</button>\n'.format(
                a=active, pid=pid, n=name
            )
        html += "</div>\n\n"

    # Go through datasets creating plots
    for pidx, pdata in enumerate(plotdata):
        # Save plot data to file
        fdata = {}
        for d in pdata:
            for didx, dval in enumerate(d["data"]):
                s_name = plotsamples[pidx][didx]
                if s_name not in fdata:
                    fdata[s_name] = dict()
                fdata[s_name][d["name"]] = dval
        if pconfig.get("save_data_file", True):
            util_functions.write_data_file(fdata, pids[pidx])

        # Plot percentage as well as counts
        plot_pcts = [False]
        if pconfig.get("cpswitch") is not False:
            plot_pcts = [False, True]

        # Switch out NaN for 0s so that MatPlotLib doesn't ignore stuff
        for idx, d in enumerate(pdata):
            pdata[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        for plot_pct in plot_pcts:
            # Plot ID
            pid = pids[pidx]
            hide_plot = False
            if plot_pct is True:
                pid = f"{pid}_pc"
                if pconfig.get("cpswitch_c_active", True) is True:
                    hide_plot = True
            else:
                if pconfig.get("cpswitch_c_active", True) is not True:
                    hide_plot = True

            # Set up figure

            # Height has a default, then adjusted by the number of samples
            plt_height = len(plotsamples[pidx]) / 2.3  # Default in inches, empirically determined
            plt_height = max(6, plt_height)  # At least 6" tall
            plt_height = min(30, plt_height)  # Cap at 30" tall

            # Use fixed height if pconfig['height'] is set (convert pixels -> inches)
            if "height" in pconfig:
                # Default interactive height in pixels = 512
                # Not perfect replication, but good enough
                plt_height = 6 * (pconfig["height"] / 512)

            bar_width = 0.8

            fig = plt.figure(figsize=(14, plt_height), frameon=False)
            axes = fig.add_subplot(111)
            y_ind = range(len(plotsamples[pidx]))

            # Count totals for each sample
            if plot_pct is True:
                s_totals = [0 for _ in pdata[0]["data"]]
                for series_idx, d in enumerate(pdata):
                    for sample_idx, v in enumerate(d["data"]):
                        s_totals[sample_idx] += v

            # Plot bars
            dlabels = []
            prev_values = None
            for idx, d in enumerate(pdata):
                # Plot percentages
                values = [x for x in d["data"]]
                if len(values) < len(y_ind):
                    values.extend([0] * (len(y_ind) - len(values)))
                if plot_pct is True:
                    for key, var in enumerate(values):
                        s_total = s_totals[key]
                        if s_total == 0:
                            values[key] = 0
                        else:
                            values[key] = (float(var + 0.0) / float(s_total)) * 100

                # Get offset for stacked bars
                if idx == 0:
                    prevdata = [0] * len(plotsamples[pidx])
                else:
                    for i, p in enumerate(prevdata):
                        prevdata[i] += prev_values[i]
                # Save the name of this series
                dlabels.append(d["name"])
                # Add the series of bars to the plot
                axes.barh(
                    y_ind,
                    values,
                    bar_width,
                    left=prevdata,
                    color=d["color"],
                    align="center",
                    linewidth=pconfig.get("borderWidth", 0),
                )
                prev_values = values

            # Tidy up axes
            axes.tick_params(
                labelsize=pconfig.get("labelSize", 8), direction="out", left=False, right=False, top=False, bottom=False
            )
            axes.set_xlabel(pconfig.get("ylab", ""))  # I know, I should fix the fact that the config is switched
            axes.set_ylabel(pconfig.get("xlab", ""))
            axes.set_yticks(y_ind)  # Specify where to put the labels
            axes.set_yticklabels(plotsamples[pidx])  # Set y axis sample name labels
            axes.set_ylim((-0.5, len(y_ind) - 0.5))  # Reduce padding around plot area
            if plot_pct is True:
                axes.set_xlim((0, 100))
                # Add percent symbols
                vals = axes.get_xticks()
                axes.set_xticks(axes.get_xticks())
                axes.set_xticklabels([f"{x:.0f}%" for x in vals])
            else:
                default_xlimits = axes.get_xlim()
                axes.set_xlim((pconfig.get("ymin", default_xlimits[0]), pconfig.get("ymax", default_xlimits[1])))
            if "title" in pconfig:
                top_gap = 1 + (0.5 / plt_height)
                plt.text(
                    0.5, top_gap, pconfig["title"], horizontalalignment="center", fontsize=16, transform=axes.transAxes
                )
            axes.grid(True, zorder=0, which="both", axis="x", linestyle="-", color="#dedede", linewidth=1)
            axes.set_axisbelow(True)
            axes.spines["right"].set_visible(False)
            axes.spines["top"].set_visible(False)
            axes.spines["bottom"].set_visible(False)
            axes.spines["left"].set_visible(False)
            plt.gca().invert_yaxis()  # y-axis is reverse sorted otherwise

            # Hide some labels if we have a lot of samples
            show_nth = max(1, math.ceil(len(pdata[0]["data"]) / 150))
            for idx, label in enumerate(axes.get_yticklabels()):
                if idx % show_nth != 0:
                    label.set_visible(False)

            # Legend
            bottom_gap = -1 * (1 - ((plt_height - 1.5) / plt_height))
            lgd = axes.legend(
                dlabels,
                loc="lower center",
                bbox_to_anchor=(0, bottom_gap, 1, 0.102),
                ncol=5,
                mode="expand",
                fontsize=pconfig.get("labelSize", 8),
                frameon=False,
            )

            # Should this plot be hidden on report load?
            hidediv = ""
            if pidx > 0 or hide_plot:
                hidediv = ' style="display:none;"'

            # Save the plot to the data directory if export is requested
            if config.export_plots:
                for fformat in config.export_plot_formats:
                    # Make the directory if it doesn't already exist
                    plot_dir = os.path.join(config.plots_dir, fformat)
                    if not os.path.exists(plot_dir):
                        os.makedirs(plot_dir)
                    # Save the plot
                    plot_fn = os.path.join(plot_dir, f"{pid}.{fformat}")
                    fig.savefig(plot_fn, format=fformat, bbox_extra_artists=(lgd,), bbox_inches="tight")

            # Output the figure to a base64 encoded string
            if getattr(get_template_mod(), "base64_plots", True) is True:
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format="png", bbox_inches="tight")
                b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
                img_buffer.close()
                html += '<div class="mqc_mplplot" id="{}"{}><img src="data:image/png;base64,{}" /></div>'.format(
                    pid, hidediv, b64_img
                )

            # Link to the saved image
            else:
                plot_relpath = os.path.join(config.plots_dir_name, "png", f"{pid}.png")
                html += f'<div class="mqc_mplplot" id="{pid}"{hidediv}><img src="{plot_relpath}" /></div>'

            plt.close(fig)

    # Close wrapping div
    html += "</div>"

    return html
