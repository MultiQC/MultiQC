#!/usr/bin/env python

""" MultiQC functions to plot a scatter plot """

import logging
import random

from multiqc.utils import config, report

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"


def plot(data, pconfig=None):
    """Plot a scatter plot with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]

    # Generate the data dict structure expected by HighCharts series
    plotdata = list()
    for data_index, ds in enumerate(data):
        d = list()
        for s_name in ds:
            # Ensure any overwritting conditionals from data_labels (e.g. ymax) are taken in consideration
            series_config = pconfig.copy()
            if (
                "data_labels" in pconfig and type(pconfig["data_labels"][data_index]) is dict
            ):  # if not a dict: only dataset name is provided
                series_config.update(pconfig["data_labels"][data_index])

            if type(ds[s_name]) is not list:
                ds[s_name] = [ds[s_name]]
            for k in ds[s_name]:
                if k["x"] is not None:
                    if "xmax" in series_config and float(k["x"]) > float(series_config["xmax"]):
                        continue
                    if "xmin" in series_config and float(k["x"]) < float(series_config["xmin"]):
                        continue
                if k["y"] is not None:
                    if "ymax" in series_config and float(k["y"]) > float(series_config["ymax"]):
                        continue
                    if "ymin" in series_config and float(k["y"]) < float(series_config["ymin"]):
                        continue
                this_series = {"x": k["x"], "y": k["y"]}
                try:
                    this_series["name"] = "{}: {}".format(s_name, k["name"])
                except KeyError:
                    this_series["name"] = s_name
                try:
                    this_series["color"] = k["color"]
                except KeyError:
                    try:
                        this_series["color"] = series_config["colors"][s_name]
                    except KeyError:
                        pass
                d.append(this_series)
        plotdata.append(d)

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

    # Make a plot
    return highcharts_scatter_plot(plotdata, pconfig)


def highcharts_scatter_plot(plotdata, pconfig=None):
    """
    Build the HTML needed for a HighCharts scatter plot. Should be
    called by scatter.plot(), which properly formats input data.
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
                xlab = 'data-xlab="{}"'.format(name) if name != k + 1 else ""
            html += '<button class="btn btn-default btn-sm {a}" data-action="set_data" {y} {ym} {xl} data-newdata="{k}" data-target="{id}">{n}</button>\n'.format(
                a=active, id=pconfig["id"], n=name, y=ylab, ym=ymax, xl=xlab, k=k
            )
        html += "</div>\n\n"

    # The plot div
    html += '<div class="hc-plot-wrapper"{height}><div id="{id}" class="hc-plot not_rendered hc-scatter-plot"><small>loading..</small></div></div></div> \n'.format(
        id=pconfig["id"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[pconfig["id"]] = {"plot_type": "scatter", "datasets": plotdata, "config": pconfig}

    return html
