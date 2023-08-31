#!/usr/bin/env python

""" MultiQC functions to plot a linegraph """

import base64
import inspect
import io
import logging
import os
import random
import re
import sys
from collections import OrderedDict

from multiqc.utils import config, mqc_colour, report, util_functions

logger = logging.getLogger(__name__)

try:
    # Import matplot lib but avoid default X environment
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logger.debug("Using matplotlib version {}".format(matplotlib.__version__))
except Exception as e:
    # MatPlotLib can break in a variety of ways. Fake an error message and continue without it if so.
    # The lack of the library will be handled when plots are attempted
    print("##### ERROR! MatPlotLib library could not be loaded!    #####", file=sys.stderr)
    print("##### Flat plots will instead be plotted as interactive #####", file=sys.stderr)
    print(e)

letters = "abcdefghijklmnopqrstuvwxyz"

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(data, pconfig=None):
    """Plot a line graph with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    # Don't just use {} as the default argument as it's mutable. See:
    # http://python-guide-pt-br.readthedocs.io/en/latest/writing/gotchas/
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]

    # Validate config if linting
    if config.lint:
        # Get module name
        modname = ""
        callstack = inspect.stack()
        for n in callstack:
            if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                callpath = n[1].split("multiqc/modules/", 1)[-1]
                modname = ">{}< ".format(callpath)
                break
        # Look for essential missing pconfig keys
        for k in ["id", "title", "ylab"]:
            if k not in pconfig:
                errmsg = "LINT: {}Linegraph pconfig was missing key '{}'".format(modname, k)
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
        # Check plot title format
        if not re.match(r"^[^:]*\S: \S[^:]*$", pconfig.get("title", "")):
            errmsg = "LINT: {} Linegraph title did not match format 'Module: Plot Name' (found '{}')".format(
                modname, pconfig.get("title", "")
            )
            logger.error(errmsg)
            report.lint_errors.append(errmsg)

    # Smooth dataset if requested in config
    if pconfig.get("smooth_points", None) is not None:
        sumcounts = pconfig.get("smooth_points_sumcounts", True)
        for i, d in enumerate(data):
            if type(sumcounts) is list:
                sumc = sumcounts[i]
            else:
                sumc = sumcounts
            data[i] = smooth_line_data(d, pconfig["smooth_points"], sumc)

    # Add sane plotting config defaults
    for idx, yp in enumerate(pconfig.get("yPlotLines", [])):
        pconfig["yPlotLines"][idx]["width"] = pconfig["yPlotLines"][idx].get("width", 2)

    # Add initial axis labels if defined in `data_labels` but not main config
    if pconfig.get("ylab") is None:
        try:
            pconfig["ylab"] = pconfig["data_labels"][0]["ylab"]
        except (KeyError, IndexError):
            pass
    if pconfig.get("xlab") is None:
        try:
            pconfig["xlab"] = pconfig["data_labels"][0]["xlab"]
        except (KeyError, IndexError):
            pass

    # Generate the data dict structure expected by HighCharts series
    plotdata = list()
    for data_index, d in enumerate(data):
        thisplotdata = list()

        for s in sorted(d.keys()):
            # Ensure any overwritten conditionals from data_labels (e.g. ymax) are taken in consideration
            series_config = pconfig.copy()
            if (
                "data_labels" in pconfig and type(pconfig["data_labels"][data_index]) is dict
            ):  # if not a dict: only dataset name is provided
                series_config.update(pconfig["data_labels"][data_index])

            pairs = list()
            maxval = 0
            if "categories" in series_config:
                if "categories" not in pconfig or type(pconfig["categories"]) is not list:
                    pconfig["categories"] = list()
                # Add any new categories
                for k in d[s].keys():
                    if k not in pconfig["categories"]:
                        pconfig["categories"].append(k)
                # Go through categories and add either data or a blank
                for k in pconfig["categories"]:
                    try:
                        pairs.append(d[s][k])
                        maxval = max(maxval, d[s][k])
                    except KeyError:
                        pairs.append(None)
            else:
                # Discard > ymax or just hide?
                # If it never comes back into the plot, discard. If it goes above then comes back, just hide.
                discard_ymax = None
                discard_ymin = None
                for k in sorted(d[s].keys()):
                    if "xmax" in series_config and float(k) > float(series_config["xmax"]):
                        continue
                    if "xmin" in series_config and float(k) < float(series_config["xmin"]):
                        continue
                    if d[s][k] is not None and "ymax" in series_config:
                        if float(d[s][k]) > float(series_config["ymax"]):
                            discard_ymax = True
                        elif discard_ymax is True:
                            discard_ymax = False
                    if d[s][k] is not None and "ymin" in series_config:
                        if float(d[s][k]) > float(series_config["ymin"]):
                            discard_ymin = True
                        elif discard_ymin is True:
                            discard_ymin = False

                # Build the plot data structure
                for k in sorted(d[s].keys()):
                    if k is not None:
                        if "xmax" in series_config and float(k) > float(series_config["xmax"]):
                            continue
                        if "xmin" in series_config and float(k) < float(series_config["xmin"]):
                            continue
                    if d[s][k] is not None:
                        if (
                            "ymax" in series_config
                            and float(d[s][k]) > float(series_config["ymax"])
                            and discard_ymax is not False
                        ):
                            continue
                        if (
                            "ymin" in series_config
                            and float(d[s][k]) < float(series_config["ymin"])
                            and discard_ymin is not False
                        ):
                            continue
                    pairs.append([k, d[s][k]])
                    try:
                        maxval = max(maxval, d[s][k])
                    except TypeError:
                        pass
            if maxval > 0 or series_config.get("hide_empty") is not True:
                this_series = {"name": s, "data": pairs}
                try:
                    this_series["color"] = series_config["colors"][s]
                except:
                    pass
                thisplotdata.append(this_series)
        plotdata.append(thisplotdata)

    # Add on annotation data series
    try:
        if pconfig.get("extra_series"):
            extra_series = pconfig["extra_series"]
            if type(pconfig["extra_series"]) == dict:
                extra_series = [[pconfig["extra_series"]]]
            elif type(pconfig["extra_series"]) == list and type(pconfig["extra_series"][0]) == dict:
                extra_series = [pconfig["extra_series"]]
            for i, es in enumerate(extra_series):
                for s in es:
                    plotdata[i].append(s)
    except (KeyError, IndexError):
        pass

    # Add colors to the categories if not set. Since the "plot_defaults" scale is
    # identical to default scale of the Highcharts JS library, this is not strictly
    # needed. But it future proofs when we replace Highcharts with something else.
    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    for si, sd in enumerate(plotdata):
        for di, d in enumerate(sd):
            d.setdefault("color", scale.get_colour(di, lighten=1))

    # Make a plot - template custom, or interactive or flat
    try:
        return get_template_mod().linegraph(plotdata, pconfig)
    except (AttributeError, TypeError):
        if config.plots_force_flat or (
            not config.plots_force_interactive and plotdata and len(plotdata[0]) > config.plots_flat_numseries
        ):
            try:
                report.num_mpl_plots += 1
                return matplotlib_linegraph(plotdata, pconfig)
            except Exception as e:
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
            except:
                name = k + 1
            try:
                ylab = 'data-ylab="{}"'.format(pconfig["data_labels"][k]["ylab"])
            except:
                ylab = 'data-ylab="{}"'.format(name) if name != k + 1 else ""
            try:
                ymax = 'data-ymax="{}"'.format(pconfig["data_labels"][k]["ymax"])
            except:
                ymax = ""
            try:
                xlab = 'data-xlab="{}"'.format(pconfig["data_labels"][k]["xlab"])
            except:
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
        except:
            name = k + 1
        pid = "mqc_{}_{}".format(pconfig["id"], name)
        pid = report.save_htmlid(pid, skiplint=True)
        pids.append(pid)

    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )
    html += '<div class="mqc_mplplot_plotgroup" id="{}">'.format(pconfig["id"])

    # Buttons to cycle through different datasets
    if len(plotdata) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(plotdata):
            pid = pids[k]
            active = "active" if k == 0 else ""
            try:
                name = pconfig["data_labels"][k]["name"]
            except:
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
            fdata = OrderedDict()
            lastcats = None
            sharedcats = True
            for d in pdata:
                fdata[d["name"]] = OrderedDict()

                # Check to see if all categories are the same
                if len(d["data"]) > 0 and type(d["data"][0]) is list:
                    if lastcats is None:
                        lastcats = [x[0] for x in d["data"]]
                    elif lastcats != [x[0] for x in d["data"]]:
                        sharedcats = False

                for i, x in enumerate(d["data"]):
                    if type(x) is list:
                        fdata[d["name"]][str(x[0])] = x[1]
                    else:
                        try:
                            fdata[d["name"]][pconfig["categories"][i]] = x
                        except (KeyError, IndexError):
                            fdata[d["name"]][str(i)] = x

            # Custom tsv output if the x-axis varies
            if not sharedcats and config.data_format == "tsv":
                fout = ""
                for d in pdata:
                    fout += "\t" + "\t".join([str(x[0]) for x in d["data"]])
                    fout += "\n{}\t".format(d["name"])
                    fout += "\t".join([str(x[1]) for x in d["data"]])
                    fout += "\n"
                with io.open(os.path.join(config.data_dir, "{}.txt".format(pid)), "w", encoding="utf-8") as f:
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

        # Tidy up axes
        axes.tick_params(
            labelsize=pconfig.get("labelSize", 8), direction="out", left=False, right=False, top=False, bottom=False
        )
        axes.set_xlabel(pconfig.get("xlab", ""))
        axes.set_ylabel(pconfig.get("ylab", ""))

        # Dataset specific y label
        try:
            axes.set_ylabel(pconfig["data_labels"][pidx]["ylab"])
        except:
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
        except:
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
            xmax = xmin + pconfig["xMinRange"]
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
                plot_fn = os.path.join(plot_dir, "{}.{}".format(pid, fformat))
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
            plot_relpath = os.path.join(config.plots_dir_name, "png", "{}.png".format(pid))
            html += '<div class="mqc_mplplot" id="{}"{}><img src="{}" /></div>'.format(pid, hidediv, plot_relpath)

        plt.close(fig)

    # Close wrapping div
    html += "</div>"

    return html


def smooth_line_data(data, numpoints, sumcounts=True):
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
    incides: [0.0, 4.5, 9] -> [0, 5, 9]
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
        smoothed_d = OrderedDict(xy for i, xy in enumerate(d.items()) if i in first_element_indices)
        smoothed_data[s_name] = smoothed_d

    return smoothed_data
