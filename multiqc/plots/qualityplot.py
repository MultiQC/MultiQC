#!/usr/bin/env python

""" Functions to plot miniature boxplots for ber pase quality """
import base64
import io
import logging
import os
import random
import sys

from multiqc.utils import config, report

logger = logging.getLogger(__name__)

try:
    # Import matplot lib but avoid default X environment
    import matplotlib
    from matplotlib.lines import Line2D
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
    """Plot per base quality distribution.
    :param data: dictionary. Keys are samples. Values are dictionaries of quality metrics per base.
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

    # Add sane plotting config defaults
    for idx, yp in enumerate(pconfig.get("yPlotLines", [])):
        pconfig["yPlotLines"][idx]["width"] = pconfig["yPlotLines"][idx].get(
            "width", 2)

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
                    data[i].append(s)
    except (KeyError, IndexError):
        pass

    # Make the plot
    return matplotlib_quality_plot(data, pconfig)


def matplotlib_quality_plot(plotdata, pconfig=None):
    if pconfig is None:
        pconfig = {}

    # Plot ID
    if pconfig.get("id") is None:
        pconfig["id"] = "mqc_mplplot_" + "".join(random.sample(letters, 10))

    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )

    # Same defaults as HighCharts for consistency
    default_colors = [
        "#7cb5ec",
        "#434348",
        "#136902",  # Darker shade of green to contrast with the light green background
        "#f7a35c",
        "#8085e9",
        "#f15c80",
        "#e4d354",
        "#2b908f",
        "#f45b5b",
        "#91e8e1",
    ]

    plt_height = 6
    # Use fixed height if pconfig['height'] is set (convert pixels -> inches)
    if "height" in pconfig:
        # Default interactive height in pixels = 512
        # Not perfect replication, but good enough
        plt_height = 6 * (pconfig["height"] / 512)

    # Set up figure
    fig = plt.figure(figsize=(16, plt_height), frameon=False)
    axes = fig.add_subplot(111)

    # Plot data
    sample_names = list(plotdata.keys())
    n_samples = len(sample_names)
    for s, i, c in zip(sample_names, range(n_samples), default_colors):
        # Solid lines : lower quartile to upper quartile
        axes.vlines([x + i/(n_samples+1) for x in list(plotdata[s].keys())],
                    [m['lower_quartile'] for m in list(plotdata[s].values())],
                    [m['upper_quartile'] for m in list(plotdata[s].values())],
                    linewidth = 1, color=c, label=s)
        # Dotted lines : 10th to 90th percentile
        axes.vlines([x + i/(n_samples+1) for x in list(plotdata[s].keys())],
                    [m['10th_percentile'] for m in list(plotdata[s].values())],
                    [m['90th_percentile'] for m in list(plotdata[s].values())],
                    linewidth = 1, linestyle='dotted', color=c)
        # Center : median
        axes.vlines([x + i/(n_samples+1) for x in list(plotdata[s].keys())],
                    [m['median'] - 0.2 for m in list(plotdata[s].values())],
                    [m['median'] + 0.2 for m in list(plotdata[s].values())],
                    linewidth = 2, linestyle='solid', color=c)

    # y axis limits
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

    # x axis limits
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
        plt.text(0.5, 1.05, pconfig["title"], horizontalalignment="center",
                 fontsize=16, transform=axes.transAxes)

    # X axis categories, if specified
    if "categories" in pconfig:
        axes.set_xticks([i for i, v in enumerate(pconfig["categories"])])
        axes.set_xticklabels(pconfig["categories"])

    # Axis lines
    xlim = axes.get_xlim()
    axes.plot([xlim[0], xlim[1]], [0, 0], linestyle="-",
              color="#dedede", linewidth=2)
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
    if len(plotdata) <= 15:
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

    # Save the plot to the data directory if export is requests
    if config.export_plots:
        for fformat in config.export_plot_formats:
            # Make the directory if it doesn't already exist
            plot_dir = os.path.join(config.plots_dir, fformat)
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
            # Save the plot
            plot_fn = os.path.join(plot_dir, "{}.{}".format(fformat))
            fig.savefig(plot_fn, format=fformat, bbox_inches="tight")

    # Output the figure to a base64 encoded string
    if getattr(get_template_mod(), "base64_plots", True) is True:
        img_buffer = io.BytesIO()
        fig.savefig(img_buffer, format="png", bbox_inches="tight")
        b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
        img_buffer.close()
        html += '<div class="mqc_mplplot" id="{}"><img src="data:image/png;base64,{}" /></div>'.format(pconfig["id"], b64_img)

    # Save to a file and link <img>
    else:
        plot_relpath = os.path.join(config.plots_dir_name, "png", "{}.png")
        html += '<div class="mqc_mplplot" id="{}"><img src="{}" /></div>'.format(pconfig["id"], plot_relpath)

    plt.close(fig)

    # Close wrapping div
    html += "</div>"

    return html
