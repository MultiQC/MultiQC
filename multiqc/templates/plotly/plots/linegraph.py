"""Plotly bargraph functionality."""
import base64
import inspect
import io
import logging
import math
import os
import re
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Union, cast

import plotly.graph_objects as go

from multiqc.templates.plotly.plots import get_template_mod
from multiqc.utils import config, mqc_colour, report, util_functions

logger = logging.getLogger(__name__)


class Settings:
    def __init__(self, pconfig: dict):
        # Validate config if linting
        if config.strict:
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

        self.title = pconfig.get("title")
        self._id = pconfig.get("id")
        self.xlab = pconfig.get("xlab")
        self.ylab = pconfig.get("ylab")

        # Add initial axis labels if defined in `data_labels` but not main config
        if self.ylab is None:
            try:
                self.ylab = pconfig["data_labels"][0]["ylab"]
            except (KeyError, IndexError):
                pass
        if self.xlab is None:
            try:
                self.xlab = pconfig["data_labels"][0]["xlab"]
            except (KeyError, IndexError):
                pass

        self.smooth_points = pconfig.get("smooth_points", None) is not None
        self.sum_counts = pconfig.get("smooth_points_sumcounts", True)
        self.y_plot_lines = pconfig.get("yPlotLines", [])
        # Add sane plotting config defaults
        for idx, yp in enumerate(self.y_plot_lines):
            self.y_plot_lines[idx]["width"] = self.y_plot_lines[idx].get("width", 2)

        self.data_labels: List = pconfig.get("data_labels", [])
        self.save_data_file = pconfig.get("save_data_file", True)

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, id: str):
        # Sanitise plot ID and check for duplicates
        self._id = report.save_htmlid(id)


def _layout(settings) -> go.Layout:
    """
    Layout object for the plot.
    """
    return go.Layout(
        title=dict(
            text=settings.title,
            xanchor="center",
            x=0.5,
            font=dict(size=20),
        ),
        yaxis=dict(
            # the plot is "transposed", so yaxis corresponds to the horizontal axis
            title=dict(text=settings.xlab),
            showgrid=False,
            categoryorder="category descending",
            automargin=True,
        ),
        xaxis=dict(
            title=dict(text=settings.ylab),
            gridcolor="rgba(0,0,0,0.1)",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode="y unified",
        # template="simple_white",
        font={"color": "Black", "family": "Lucida Grande"},
        colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
        barmode="stack",
        autosize=True,
        margin=dict(
            pad=10,  # pad sample names a bit
        ),
        annotations=[
            dict(
                text="Created with MultiQC",
                font=dict(size=10, color="rgba(0,0,0,0.5)"),
                xanchor="right",
                yanchor="top",
                x=1.05,
                y=1.05,
                textangle=-90,
                # yanchor="bottom",
                # x=1.07,
                # y=-0.35,
                xref="paper",
                yref="paper",
                showarrow=False,
            )
        ],
        modebar=dict(
            bgcolor="rgba(0, 0, 0, 0)",
            color="rgba(0, 0, 0, 0.5)",
            activecolor="rgba(0, 0, 0, 1)",
        ),
    )


def _figure(
    layout: go.Layout,
    data: List[Dict],
) -> go.Figure:
    """
    Build one Plotly Figure object.
    """
    fig = go.Figure()
    fig.update_layout(
        # The function expects a dict, even though go.Layout works just fine
        cast(dict, layout),
    )
    # for data in data_by_cat:
    #     fig.add_trace(
    #         go.Bar(
    #             y=sample_names,
    #             x=data["data"],
    #             name=data["name"],
    #             orientation="h",
    #             marker=dict(color=data.get("color"), line=dict(width=0)),
    #         ),
    #     )
    return fig


def linegraph(
    data: Union[List[Dict], Dict],
    pconfig: Dict,
) -> str:
    return add_to_report(data, Settings(pconfig))


def add_to_report(
    data: Union[List[Dict], Dict],
    settings: Settings,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    # Given one dataset - turn it into a list
    if type(data) is not list:
        data = [data]

    # Smooth dataset if requested in config
    if settings.smooth_points:
        for i, d in enumerate(data):
            data[i] = _smooth_line_data(d, settings.smooth_points)

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
    mod = get_template_mod()
    if "linegraph" in mod.__dict__ and callable(mod.linegraph):
        try:
            return mod.linegraph(plotdata, pconfig)
        except:
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise
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

    if config.plots_force_flat:
        return _add_to_report_flat(
            settings,
            datas,
        )
    else:
        return _add_to_report_interactive(
            settings,
            datas,
        )


def _smooth_line_data(data, numpoints):
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
        smoothed_d = {xy for i, xy in enumerate(d.items()) if i in first_element_indices}
        smoothed_data[s_name] = smoothed_d

    return smoothed_data


def _add_to_report_interactive(
    settings: Settings,
    data_by_cat_lists: List[List[Dict]],
    samples_lists: List[List[str]],
) -> str:
    html = '<div class="mqc_hcplot_plotgroup">'
    btn_tmpl = (
        f'<button class="{{cls}} btn btn-default btn-sm {{active}}" data-target="{self.id}">{{label}}' f"</button> \n"
    )
    # Counts / Percentages / Log Switches
    if settings.add_pct_tab or settings.add_log_tab:
        # html += '<div class="btn-group hc_switch_group"> \n'
        if settings.add_pct_tab:
            html += btn_tmpl.format(
                active=settings.p_active,
                label=settings.p_label,
                cls="switch_percent",
            )
        if settings.add_log_tab:
            html += btn_tmpl.format(
                active=settings.l_active,
                label=settings.l_label,
                cls="switch_log10",
            )
        # html += "</div> "
        if len(data_by_cat_lists) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(data_by_cat_lists) > 1:
        html += '<div class="btn-group dataset_switch_group">\n'
        for k, p in enumerate(data_by_cat_lists):
            active = "active" if k == 0 else ""
            try:
                name = settings.data_labels[k]["name"]
            except (TypeError, KeyError):
                try:
                    name = settings.data_labels[k]
                except KeyError:
                    name = k + 1
            try:
                ymax = f'data-ymax="{settings.data_labels[k]["ymax"]}"'
            except:
                ymax = ""
            html += f'<button class="btn btn-default btn-sm {active}" {ymax} data-dataset_index="{k}" data-target="{self.id}">{name}</button>\n'
        html += "</div>\n\n"

    # Plot HTML
    html += """<div class="hc-plot-wrapper"{height}>
        <div id="{id}" class="hc-plot not_rendered hc-bar-plot"></div>
    </div></div>""".format(
        id=settings.id,
        height=f' style="height:{settings.height}px"' if settings.height else "",
    )

    report.num_hc_plots += 1

    report.plot_data[settings.id] = {
        "plot_type": "bar_graph",
        "samples": samples_lists,
        "datasets": data_by_cat_lists,
        "layout": _layout(settings).to_plotly_json(),
        "active_dataset_idx": 0,
        "p_active": settings.p_active,
        "l_active": settings.l_active,
    }

    return html


def _add_to_report_flat(
    settings: Settings,
    datas: List[Dict],
) -> str:
    """
    Create an HTML for a static plot version.
    """
    # Individual plot IDs
    pids = []
    for k in range(len(datas)):
        try:
            name = settings.data_labels[k]
        except IndexError:
            name = k + 1
        pid = report.save_htmlid(f"{settings.id}_{name}", skiplint=True)
        pids.append(pid)

    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )
    html += f'<div class="mqc_mplplot_plotgroup" id="{settings.id}">'

    # Counts / Percentages Switch
    if settings.add_pct_tab and not config.simple_output:
        html += f'<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_setcountspcnt"> \n\
            <button class="btn btn-default btn-sm {settings.c_active} counts">{settings.c_label}</button> \n\
            <button class="btn btn-default btn-sm {settings.p_active} pcnt">{settings.p_label}</button> \n\
        </div> '
        if len(data_by_cat_lists) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(data_by_cat_lists) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(data_by_cat_lists):
            pid = pids[k]
            active = "active" if k == 0 else ""
            try:
                name = settings.data_labels[k]
            except:
                name = k + 1
            html += f'<button class="btn btn-default btn-sm {active}" data-target="#{pid}">{name}</button>\n'
        html += "</div>\n\n"

        # Collect all categories, and fill in with zeroes for samples that having
        # any of cats missing
        cat_to_color = dict()
        for p in data_by_cat_lists:
            for d in p:
                cat_to_color[d["name"]] = d["color"]
        for data_by_cat in data_by_cat_lists:
            for cat, color in cat_to_color.items():
                if cat not in [d["name"] for d in data_by_cat]:
                    data_by_cat.append(
                        {
                            "name": cat,
                            "color": color,
                            "data": [0 for _ in samples_lists[0]],
                            "data_pct": [0.0 for _ in samples_lists[0]],
                            "data_log": [0.0 for _ in samples_lists[0]],
                        }
                    )
        # Sort categories by name
        for data_by_cat in data_by_cat_lists:
            data_by_cat.sort(key=lambda x: x["name"])

    # Finally, build and save plots
    for pidx, data in enumerate(datas):
        # Save plot data to file
        fdata = {}
        # for d in data_by_cat:
        #     for didx, dval in enumerate(d["data"]):
        #         s_name = samples[didx]
        #         if s_name not in fdata:
        #             fdata[s_name] = dict()
        #         fdata[s_name][d["name"]] = dval
        if settings.save_data_file:
            util_functions.write_data_file(fdata, pids[pidx])

        # Switch out NaN for 0s
        for idx, d in enumerate(data_by_cat):
            data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        # Calculate log10 values
        if settings.add_log_tab:
            for cat_idx, d in enumerate(data_by_cat):
                values = [x for x in d["data"]]
                # if len(values) < len(samples):
                #     values.extend([0] * (len(samples) - len(values)))
                for key, var in enumerate(values):
                    if var == 0:
                        values[key] = 0
                    else:
                        values[key] = math.log10(var)
                d["data_log"] = values

        View = namedtuple(
            "View",
            [
                "data",
                "active",
                "suffix",
                "label",
                "xaxis_tickformat",
            ],
        )
        views = [
            View(
                data,
                active=not settings.p_active,
                suffix="",
                label=settings.c_label,
                xaxis_tickformat="",
            ),
        ]
        # if settings.add_pct_tab:
        #     views.append(
        #         View(
        #             [
        #                 {
        #                     "data": data_by_cat[i]["data_pct"],
        #                     "name": data_by_cat[i]["name"],
        #                     "color": data_by_cat[i]["color"],
        #                 }
        #                 for i, vals in enumerate(data_by_cat)
        #             ],
        #             active=settings.p_active,
        #             suffix="_pc",
        #             label=settings.p_label,
        #             xaxis_tickformat=".0%",
        #         )
        #     )
        # if settings.add_log_tab:
        #     views.append(
        #         View(
        #             [
        #                 {
        #                     "data": data_by_cat[i]["data_log"],
        #                     "name": data_by_cat[i]["name"],
        #                     "color": data_by_cat[i]["color"],
        #                 }
        #                 for i, vals in enumerate(data_by_cat)
        #             ],
        #             active=False,
        #             suffix="_log",
        #             label=settings.l_label,
        #             xaxis_tickformat="",
        #         )
        #     )

        for view in views:
            fig = _figure(
                _layout(settings),
                data=view.data,
            )
            if view.xaxis_tickformat:
                fig.update_layout({"xaxis": {"tickformat": view.xaxis_tickformat}})

            # Should this plot be hidden on report load?
            hide_div = ""
            if pidx > 0 or not view.active:
                hide_div = ' style="display:none;"'

            pid = f"{pids[pidx]}{view.suffix}"

            # Save the plot to the data directory if export is requested
            if config.export_plots:
                for fformat in config.export_plot_formats:
                    # Make the directory if it doesn't already exist
                    plot_dir = os.path.join(config.plots_dir, fformat)
                    if not os.path.exists(plot_dir):
                        os.makedirs(plot_dir)
                    # Save the plot
                    plot_fn = os.path.join(plot_dir, f"{pid}.{fformat}")
                    fig.write_image(
                        plot_fn,
                        format=fformat,
                        width=fig.layout.width,
                        height=fig.layout.height,
                        scale=1,
                    )

            # Output the figure to a base64 encoded string
            if getattr(get_template_mod(), "base64_plots", True) is True:
                img_buffer = io.BytesIO()
                fig.write_image(
                    img_buffer,
                    format="png",
                    width=fig.layout.width,
                    height=fig.layout.height,
                    scale=1,
                )
                b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
                img_buffer.close()
                html += (
                    f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="data:image/png;base64,{b64_img}" /></div>'
                )

            # Link to the saved image
            else:
                plot_relpath = Path(config.plots_dir_name) / "png" / f"{pid}.png"
                plot_relpath.parent.mkdir(parents=True, exist_ok=True)
                fig.write_image(
                    plot_relpath,
                    format="png",
                    width=900,
                    height=fig.layout.height,
                    scale=1,
                )
                html += f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'

    # Close wrapping div
    html += "</div>"
    return html
