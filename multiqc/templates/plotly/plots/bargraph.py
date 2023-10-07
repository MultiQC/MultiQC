"""Plotly bargraph functionality."""
import base64
import io
import math
import os
import random
import string
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, cast

import plotly.graph_objects as go

from multiqc.templates.plotly.plots import PlotSettings, basic_figure, get_template_mod
from multiqc.utils import config, mqc_colour, report, util_functions

"""
Currently, we have to implement the plots twice: interactive ones with highcharts 
and static versions with matplotlib.

With Plotly, we will be able to build a plot once, and then export it to both
interactive and static versions with fig.to_html() (or fig.to_plotly_json() followed
by JavaScript rendering), and with fig.write_image().

However, for static images, we would still need some interactivity, namely the buttons
to switch between percentages/logarithmic scales and counts, as well as support plots
for multiple datasets. With matplotlib, this is achieved by creating a separate
figure for each dataset, and then using the `subplots` functionality to combine them
"""


def barplot_layout(settings: PlotSettings) -> go.Layout:
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
        # height=settings.height,
        autosize=True,
        annotations=[
            dict(
                text="Created with MultiQC",
                font=dict(size=10, color="rgba(0,0,0,0.5)"),
                xanchor="right",
                yanchor="bottom",
                x=1.07,
                y=-0.35,
                xref="paper",
                yref="paper",
                showarrow=False,
            )
        ],
    )


def _bargraph_plotly_plot(
    layout: go.Layout,
    data_by_cat: List[Dict],
    sample_names: List[str],
) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        # The function expects a dict, even though go.Layout works just fine
        cast(dict, layout),
    )
    for data in data_by_cat:
        fig.add_trace(
            go.Bar(
                y=sample_names,
                x=data["data"],
                name=data["name"],
                orientation="h",
                marker=dict(color=data.get("color"), line=dict(width=0)),
            ),
        )
    return fig


def plotly_bargraph(
    data_by_cat_lists: List[List[Dict]],
    samples_lists: List[List[str]],
    pconfig: Dict,
) -> str:
    """
    :param data_by_cat_lists: List of lists of dicts with the keys: {name, color, data},
        where `name` is the category name, `color` is the color of the bar,
        and `data` is a list of values for each sample. Each outer list will
        correspond a separate tab.
    :param samples_lists: List of lists of bar names (i.e., sample names). Similarly,
        each outer list will correspond to a separate tab.
    :param pconfig: Plot parameters.
    :return: Plotly HTML
    """
    settings = PlotSettings(pconfig)

    if not settings.height:
        max_n_samples = max(len(samples) for samples in samples_lists)
        # Height has a default, then adjusted by the number of samples
        settings.height = max_n_samples // 186  # Default, empirically determined
        settings.height = max(700, settings.height)
        settings.height = min(2560, settings.height)

    uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
    is_static_suf = "static_" if config.plots_force_flat else ""
    if settings.id is None:  # ID of the plot group
        settings.id = f"mqc_{is_static_suf}plot_{uniq_suffix}"

    for pidx, (samples, data_by_cat) in enumerate(zip(samples_lists, data_by_cat_lists)):
        # Switch out NaN for 0s
        for idx, d in enumerate(data_by_cat):
            data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        # Calculate and save percentages
        if settings.add_pct_tab:
            # Count totals for each category
            sums = [0 for _ in data_by_cat[0]["data"]]
            for cat_idx, d in enumerate(data_by_cat):
                for sample_idx, v in enumerate(d["data"]):
                    if not math.isnan(v):
                        sums[sample_idx] += v
            # Now, calculate percentages for each category
            for cat_idx, d in enumerate(data_by_cat):
                values = [x for x in d["data"]]
                if len(values) < len(samples):
                    values.extend([0] * (len(samples) - len(values)))
                for key, var in enumerate(values):
                    sum_for_cat = sums[key]
                    if sum_for_cat == 0:
                        values[key] = 0
                    else:
                        values[key] = (float(var + 0.0) / float(sum_for_cat)) * 100
                d["data_pct"] = values

    if config.plots_force_flat:
        return _static_bargraph(
            data_by_cat_lists,
            samples_lists,
            settings,
        )

    html = '<div class="mqc_hcplot_plotgroup">'

    btn_tmpl = (
        f'<button class="{{cls}} btn btn-default btn-sm {{active}}" data-target="{settings.id}">{{label}}'
        f"</button> \n"
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
            html += f'<button class="btn btn-default btn-sm {active}" {ymax} data-dataset_index="{k}" data-target="{settings.id}">{name}</button>\n'
        html += "</div>\n\n"

    # Plot HTML
    html += """<div class="container-fluid hc-plot-wrapper"{height}>
        <div id="{id}" class="hc-plot not_rendered hc-bar-plot"><small>loading..</small></div>
    </div></div>""".format(
        id=settings.id,
        height=f' style="height:{settings.height}px"' if settings.height else "",
    )

    report.num_hc_plots += 1

    report.plot_data[settings.id] = {
        "plot_type": "bar_graph",
        "samples": samples_lists,
        "datasets": data_by_cat_lists,
        "config": pconfig,
        "layout": barplot_layout(settings).to_plotly_json(),
        "active_dataset_idx": 0,
        "p_active": settings.p_active,
        "l_active": settings.l_active,
    }

    return html


def _static_bargraph(
    data_by_cat_lists: List,
    samples_lists: List,
    settings: PlotSettings,
):
    # Individual plot IDs
    pids = []
    for k in range(len(data_by_cat_lists)):
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
    for pidx, (samples, data_by_cat) in enumerate(zip(samples_lists, data_by_cat_lists)):
        # Save plot data to file
        fdata = {}
        for d in data_by_cat:
            for didx, dval in enumerate(d["data"]):
                s_name = samples[didx]
                if s_name not in fdata:
                    fdata[s_name] = dict()
                fdata[s_name][d["name"]] = dval
        if settings.save_data_file:
            util_functions.write_data_file(fdata, pids[pidx])

        # Switch out NaN for 0s
        for idx, d in enumerate(data_by_cat):
            data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        # Calculate log10 values
        if settings.add_log_tab:
            for cat_idx, d in enumerate(data_by_cat):
                values = [x for x in d["data"]]
                if len(values) < len(samples):
                    values.extend([0] * (len(samples) - len(values)))
                for key, var in enumerate(values):
                    if var == 0:
                        values[key] = 0
                    else:
                        values[key] = math.log10(var)
                d["data_log"] = values

        View = namedtuple(
            "View",
            [
                "values_by_cat",
                "active",
                "suffix",
                "label",
            ],
        )
        views = [
            View(
                data_by_cat,
                active=not settings.p_active,
                suffix="",
                label=settings.c_label,
            ),
        ]
        if settings.add_pct_tab:
            views.append(
                View(
                    [
                        {
                            "data": data_by_cat[i]["data_pct"],
                            "name": data_by_cat[i]["name"],
                            "color": data_by_cat[i]["color"],
                        }
                        for i, vals in enumerate(data_by_cat)
                    ],
                    active=settings.p_active,
                    suffix="_pc",
                    label=settings.p_label,
                )
            )
        if settings.add_log_tab:
            views.append(
                View(
                    [
                        {
                            "data": data_by_cat[i]["data_log"],
                            "name": data_by_cat[i]["name"],
                            "color": data_by_cat[i]["color"],
                        }
                        for i, vals in enumerate(data_by_cat)
                    ],
                    active=False,
                    suffix="_log",
                    label=settings.l_label,
                )
            )

        for view in views:
            plot = _bargraph_plotly_plot(
                barplot_layout(settings),
                data_by_cat=view.values_by_cat,
                sample_names=samples_lists[pidx],
            )

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
                    plot.write_image(
                        plot_fn,
                        format=fformat,
                        width=plot.layout.width,
                        height=plot.layout.height,
                        scale=1,
                    )

            # Output the figure to a base64 encoded string
            if getattr(get_template_mod(), "base64_plots", True) is True:
                img_buffer = io.BytesIO()
                plot.write_image(
                    img_buffer,
                    format="png",
                    width=plot.layout.width,
                    height=plot.layout.height,
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
                plot.write_image(
                    plot_relpath,
                    format="png",
                    width=900,
                    height=plot.layout.height,
                    scale=1,
                )
                html += f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'

    # Close wrapping div
    html += "</div>"
    return html
