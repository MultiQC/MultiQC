#!/usr/bin/env python

""" MultiQC functions to plot a scatter plot """

import base64
import io
import logging
import os
import random

from multiqc.utils import config, report

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
    """Plot a box-and-whisker plot
    :param data: 2D dict, first keys as read positions, then as quantile:QV pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    if pconfig is None:
        pconfig = {}

    # Make a plot
    return matplotlib_boxplot(data, pconfig)


def mock_data(data):
    res = []
    value = 2  # Default to our minimum/ambiguous base QV
    for key in [1, 2, 5, 10, 25, 50, 75, 90]:
        value = data.get(key, value)
        if key in [10, 90]:
            res += [value] * 10
        elif key in [25, 75]:
            res += [value] * 25
        elif key == 50:
            res += [value] * 30
    return res


def mock_dataset(data):
    res = []
    for key in sorted(data.keys()):
        res.append(mock_data(data[key]))
    return res


def matplotlib_boxplot(plotdata, pconfig=None):
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
    for pidx, (pname, pdata) in enumerate(plotdata.items()):
        # Plot ID
        pid = pids[pidx]

        # Set up figure
        fig = plt.figure(figsize=(14, 6), frameon=False)
        axes = fig.add_subplot(111)
        plt.xticks(rotation=90)

        # Go through data series
        n_boxes = len(pdata)
        mock_ds = mock_dataset(pdata)
        box = axes.boxplot(mock_ds, whis=[0, 100], patch_artist=True)

        for patch in box["boxes"]:
            patch.set_facecolor("yellow")

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
            if len(plotdata) > 1:
                title = "{} for {}".format(pconfig["title"], pname)
            else:
                title = pconfig["title"]
            plt.text(0.5, 1.05, title, horizontalalignment="center", fontsize=16, transform=axes.transAxes)
        axes.set_xlabel(pconfig.get("xlab", ""))
        axes.set_ylabel(pconfig.get("ylab", ""))
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
                fontsize=8,
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

    report.num_mpl_plots += 1

    return html
