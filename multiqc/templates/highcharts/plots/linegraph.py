#!/usr/bin/env python

""" MultiQC functions to plot a linegraph """

import base64
import io
import logging
import os
import random
import sys
from typing import Dict

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


def plot(plotdata, pconfig=None):
    if config.plots_force_flat or (
        not config.plots_force_interactive and plotdata and len(plotdata[0]) > config.plots_flat_numseries
    ):
        try:
            report.num_mpl_plots += 1
            return matplotlib_linegraph(plotdata, pconfig)
        except Exception as e:
            if config.strict:
                raise
            logger.error("############### Error making MatPlotLib figure! Falling back to HighCharts.")
            logger.debug(e, exc_info=True)
            return highcharts_linegraph(plotdata, pconfig)
    else:
        # Use MatPlotLib to generate static plots if requested
        if config.export_plots:
            matplotlib_linegraph(plotdata, pconfig)
        # Return HTML for HighCharts dynamic plot
        return highcharts_linegraph(plotdata, pconfig)


def highcharts_linegraph(plotdata, pconfig=None):
    """
    Build the HTML needed for a HighCharts line graph. Should be
    called by linegraph.plot(), which properly formats input data.
    """
    if pconfig is None:
        pconfig = {}

    # Get the plot ID
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_hcplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    # Build the HTML for the page
    html = '<div class="mqc_hcplot_plotgroup">'

    # Log Switch
    if pconfig.get("logswitch") is True:
        c_active = "active"
        l_active = ""
        if pconfig.get("logswitch_active") is True:
            c_active = ""
            l_active = "active"
        c_label = pconfig.get("cpswitch_counts_label", "Counts")
        l_label = pconfig.get("logswitch_label", "Log10")
        html += '<div class="btn-group hc_switch_group"> \n'
        html += '<button class="btn btn-default btn-sm {c_a}" data-action="set_numbers" data-target="{id}" data-ylab="{c_l}">{c_l}</button> \n'.format(
            id=pconfig["id"], c_a=c_active, c_l=c_label
        )
        if pconfig.get("logswitch") is True:
            html += '<button class="btn btn-default btn-sm {l_a}" data-action="set_log" data-target="{id}" data-ylab="{l_l}">{l_l}</button> \n'.format(
                id=pconfig["id"], l_a=l_active, l_l=l_label
            )
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
            try:
                xlab = f"data-xlab=\"{pconfig['data_labels'][k]['xlab']}\""
            except Exception:
                xlab = ""
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} {x} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(
                a=active, id=pconfig["id"], n=name, y=ylab, ym=ymax, x=xlab, k=k
            )
        html += "</div>\n\n"

    # The plot div
    html += '<div class="hc-plot-wrapper"{height}><div id="{id}" class="hc-plot not_rendered hc-line-plot"><small>loading..</small></div></div></div> \n'.format(
        id=pconfig["id"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[pconfig["id"]] = {"plot_type": "xy_line", "datasets": plotdata, "config": pconfig}

    return html


def matplotlib_linegraph(plotdata, pconfig=None):
    """
    Plot a line graph with Matplot lib and return a HTML string. Either embeds a base64
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
            name = pconfig["data_labels"][k]["name"]
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

    # Buttons to cycle through different datasets
    if len(plotdata) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(plotdata):
            pid = pids[k]
            active = "active" if k == 0 else ""
            try:
                name = pconfig["data_labels"][k]["name"]
            except Exception:
                name = k + 1
            html += '<button class="btn btn-default btn-sm {a}" data-target="#{pid}">{n}</button>\n'.format(
                a=active, pid=pid, n=name
            )
        html += "</div>\n\n"

    # Go through datasets creating plots
    for pidx, pdata in enumerate(plotdata):
        # Plot ID
        pid = pids[pidx]

        # Save plot data to file
        if pconfig.get("save_data_file", True):
            fdata = dict()
            lastcats = None
            sharedcats = True
            for d in pdata:
                fdata[d["name"]] = dict()

                # Check to see if all categories are the same
                if len(d["data"]) > 0 and isinstance(d["data"][0], list):
                    if lastcats is None:
                        lastcats = [x[0] for x in d["data"]]
                    elif lastcats != [x[0] for x in d["data"]]:
                        sharedcats = False

                for i, x in enumerate(d["data"]):
                    if isinstance(x, list):
                        fdata[d["name"]][x[0]] = x[1]
                    else:
                        try:
                            fdata[d["name"]][pconfig["categories"][i]] = x
                        except Exception:
                            fdata[d["name"]][str(i)] = x

            # Custom tsv output if the x-axis varies
            if not sharedcats and config.data_format == "tsv":
                fout = ""
                for d in pdata:
                    fout += "\t" + "\t".join([str(x[0]) for x in d["data"]])
                    fout += f"\n{d['name']}\t"
                    fout += "\t".join([str(x[1]) for x in d["data"]])
                    fout += "\n"
                with io.open(os.path.join(config.data_dir, f"{pid}.txt"), "w", encoding="utf-8") as f:
                    print(fout.encode("utf-8", "ignore").decode("utf-8"), file=f)
            else:
                util_functions.write_data_file(fdata, pid)

        plt_height = 6
        # Use fixed height if pconfig['height'] is set (convert pixels -> inches)
        if "height" in pconfig:
            # Default interactive height in pixels = 512
            # Not perfect replication, but good enough
            plt_height = 6 * (pconfig["height"] / 512)

        # Set up figure
        fig = plt.figure(figsize=(14, plt_height), frameon=False)
        axes = fig.add_subplot(111)

        # Go through data series
        for idx, d in enumerate(pdata):
            # Line style
            linestyle = "solid"
            if d.get("dashStyle", None) == "Dash":
                linestyle = "dashed"

            # Reformat data (again)
            try:
                axes.plot(
                    [x[0] for x in d["data"]],
                    [x[1] for x in d["data"]],
                    label=d["name"],
                    color=d["color"],
                    linestyle=linestyle,
                    linewidth=1,
                    marker=None,
                )
            except TypeError:
                # Categorical data on x axis
                axes.plot(d["data"], label=d["name"], color=d["color"], linewidth=1, marker=None)

        # Log scale
        if pconfig.get("xLog", False):
            axes.set_xscale("log")
        if pconfig.get("yLog", False):
            axes.set_yscale("log")

        # Tidy up axes
        axes.tick_params(
            labelsize=pconfig.get("labelSize", 8), direction="out", left=False, right=False, top=False, bottom=False
        )
        axes.set_xlabel(pconfig.get("xlab", ""))
        axes.set_ylabel(pconfig.get("ylab", ""))

        # Dataset specific y label
        try:
            axes.set_ylabel(pconfig["data_labels"][pidx]["ylab"])
        except Exception:
            pass

        # Axis limits
        default_ylimits = axes.get_ylim()
        ymin = default_ylimits[0]
        if "ymin" in pconfig:
            ymin = pconfig["ymin"]
        elif "yFloor" in pconfig:
            ymin = max(pconfig["yFloor"], default_ylimits[0])
        ymax = default_ylimits[1]
        if "ymax" in pconfig:
            ymax = pconfig["ymax"]
        elif "yCeiling" in pconfig:
            ymax = min(pconfig["yCeiling"], default_ylimits[1])
        if (ymax - ymin) < pconfig.get("yMinRange", 0):
            ymax = ymin + pconfig["yMinRange"]
        axes.set_ylim((ymin, ymax))

        # Dataset specific ymax
        try:
            axes.set_ylim((ymin, pconfig["data_labels"][pidx]["ymax"]))
        except Exception:
            pass

        default_xlimits = axes.get_xlim()
        xmin = default_xlimits[0]
        if "xmin" in pconfig:
            xmin = pconfig["xmin"]
        elif "xFloor" in pconfig:
            xmin = max(pconfig["xFloor"], default_xlimits[0])
        xmax = default_xlimits[1]
        if "xmax" in pconfig:
            xmax = pconfig["xmax"]
        elif "xCeiling" in pconfig:
            xmax = min(pconfig["xCeiling"], default_xlimits[1])
        if (xmax - xmin) < pconfig.get("xMinRange", 0):
            xmax = xmin + pconfig.get("xMinRange", 0)
        axes.set_xlim((xmin, xmax))

        # Plot title
        if "title" in pconfig:
            plt.text(0.5, 1.05, pconfig["title"], horizontalalignment="center", fontsize=16, transform=axes.transAxes)
        axes.grid(True, zorder=10, which="both", axis="y", linestyle="-", color="#dedede", linewidth=1)

        # X axis categories, if specified
        if "categories" in pconfig:
            axes.set_xticks([i for i, v in enumerate(pconfig["categories"])])
            axes.set_xticklabels(pconfig["categories"])

        # Axis lines
        xlim = axes.get_xlim()
        axes.plot([xlim[0], xlim[1]], [0, 0], linestyle="-", color="#dedede", linewidth=2)
        axes.set_axisbelow(True)
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.spines["bottom"].set_visible(False)
        axes.spines["left"].set_visible(False)

        # Background colours, if specified
        if "yPlotBands" in pconfig:
            xlim = axes.get_xlim()
            for pb in pconfig["yPlotBands"]:
                axes.barh(
                    pb["from"],
                    xlim[1],
                    height=pb["to"] - pb["from"],
                    left=xlim[0],
                    color=pb["color"],
                    linewidth=0,
                    zorder=0,
                    align="edge",
                )
        if "xPlotBands" in pconfig:
            ylim = axes.get_ylim()
            for pb in pconfig["xPlotBands"]:
                axes.bar(
                    pb["from"],
                    ylim[1],
                    width=pb["to"] - pb["from"],
                    bottom=ylim[0],
                    color=pb["color"],
                    linewidth=0,
                    zorder=0,
                    align="edge",
                )

        # Tight layout - makes sure that legend fits in and stuff
        if len(pdata) <= 15:
            axes.legend(
                loc="lower center",
                bbox_to_anchor=(0, -0.22, 1, 0.102),
                ncol=5,
                mode="expand",
                fontsize=pconfig.get("labelSize", 8),
                frameon=False,
            )
            plt.tight_layout(rect=[0, 0.08, 1, 0.92])
        else:
            plt.tight_layout(rect=[0, 0, 1, 0.92])

        # Should this plot be hidden on report load?
        hidediv = ""
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
                plot_fn = os.path.join(plot_dir, f"{pid}.{fformat}")
                fig.savefig(plot_fn, format=fformat, bbox_inches="tight")

        # Output the figure to a base64 encoded string
        if getattr(get_template_mod(), "base64_plots", True) is True:
            img_buffer = io.BytesIO()
            fig.savefig(img_buffer, format="png", bbox_inches="tight")
            b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
            img_buffer.close()
            html += '<div class="mqc_mplplot" id="{}"{}><img src="data:image/png;base64,{}" /></div>'.format(
                pid, hidediv, b64_img
            )

        # Save to a file and link <img>
        else:
            plot_relpath = os.path.join(config.plots_dir_name, "png", f"{pid}.png")
            html += f'<div class="mqc_mplplot" id="{pid}"{hidediv}><img src="{plot_relpath}" /></div>'

        plt.close(fig)

    # Close wrapping div
    html += "</div>"

    return html


def smooth_line_data(data: Dict[str, Dict], numpoints: int) -> Dict[str, Dict[int, int]]:
    """
    Function to take an x-y dataset and use binning to smooth to a maximum number of datapoints.
    Each datapoint in a smoothed dataset corresponds to the first point in a bin.

    Examples to show the idea:

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=6
    we want to keep the first and the last element, thus excluding the last element from the binning:
    binsize = len([0 1 2 3 4 5 6 7 8]))/(numpoints-1) = 9/5 = 1.8
    taking points in indices rounded from multiples of 1.8: [0, 1.8, 3.6, 5.4, 7.2, 9],
    ...which evaluates to first_element_in_bin_indices=[0, 2, 4, 5, 7, 9]
    picking up the elements: [0 _ 2 _ 4 5 _ 7 _ 9]

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=9
    binsize = 9/8 = 1.125
    indices: [0.0, 1.125, 2.25, 3.375, 4.5, 5.625, 6.75, 7.875, 9] -> [0, 1, 2, 3, 5, 6, 7, 8, 9]
    picking up the elements: [0 1 2 3 _ 5 6 7 8 9]

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=3
    binsize = len(d)/numpoints = 9/2 = 4.5
    indices: [0.0, 4.5, 9] -> [0, 5, 9]
    picking up the elements: [0 _ _ _ _ 5 _ _ _ 9]
    """
    smoothed_data = dict()
    for s_name, d in data.items():
        # Check that we need to smooth this data
        if len(d) <= numpoints or len(d) == 0:
            smoothed_data[s_name] = d
            continue

        binsize = (len(d) - 1) / (numpoints - 1)
        first_element_indices = [round(binsize * i) for i in range(numpoints)]
        smoothed_d = {x: y for i, (x, y) in enumerate(d.items()) if i in first_element_indices}
        smoothed_data[s_name] = smoothed_d

    return smoothed_data
