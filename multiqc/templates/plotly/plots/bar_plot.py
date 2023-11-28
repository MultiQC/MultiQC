"""Plotly bargraph functionality."""
import base64
import io
import logging
import os
import random
import string
from pathlib import Path
from typing import Dict, List

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots import get_template_mod
from multiqc.templates.plotly.plots.plot import View, PConfig, base_layout
from multiqc.utils import config, util_functions

logger = logging.getLogger(__name__)


def bar_plot(
    datasets: List[List[Dict]],
    samples_lists: List[List[str]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a list of dicts with the keys: {name, color, data},
        where `name` is the category name, `color` is the color of the bar,
        and `data` is a list of values for each sample. Each outer list will
        correspond a separate tab.
    :param samples_lists: list of lists of bar names (that is, sample names). Similarly,
        each outer list will correspond to a separate tab.
    :param pconfig: Plot configuration dictionary
    :return: Plotly HTML
    """
    from multiqc.utils import report

    return add_to_report(
        datasets=datasets,
        samples_lists=samples_lists,
        pconfig=PConfig(pconfig, len(datasets)),
        report=report,
    )


def add_to_report(
    datasets: List[List[Dict]],
    samples_lists: List[List[str]],
    pconfig: PConfig,
    report,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    if not pconfig.height:
        max_n_samples = max(len(samples) for samples in samples_lists)
        # Height has a default, then adjusted by the number of samples
        pconfig.height = max_n_samples // 186  # Default, empirically determined
        pconfig.height = max(600, pconfig.height)
        pconfig.height = min(2560, pconfig.height)

    uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
    is_static_suf = "static_" if config.plots_force_flat else ""
    pid = report.save_htmlid(pconfig.id)
    if pid is None:  # ID of the plot group
        pid = report.save_htmlid(f"mqc_{is_static_suf}plot_{uniq_suffix}")

    # Calculate and save percentages
    if pconfig.add_pct_tab:
        for pidx, (samples, data_by_cat) in enumerate(zip(samples_lists, datasets)):
            # # Switch out NaN for 0s
            # for idx, d in enumerate(data_by_cat):
            #     data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]
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
                        values[key] = float(var + 0.0) / float(sum_for_cat)
                d["data_pct"] = values

    if config.plots_force_flat:
        html = _datasets_to_flat_imgs(
            datasets,
            samples_lists,
            pconfig,
            pid,
            # Can't use interactivity, so we will have to generate separate flat images for each
            # dataset and view. So have to make sure the individual image IDs are unique across the report:
            pids=[report.save_htmlid(f"{pid}_{dl['name']}", skiplint=True) for dl in pconfig.data_labels],
        )
    else:
        html = _datasets_to_interactive_imgs(
            datasets,
            pconfig,
            pid,
        )

        # Saving compressed data for JavaScript to pick up and uncompress.
        report.num_hc_plots += 1
        report.plot_data[pid] = {
            "plot_type": "bar_graph",
            "samples": samples_lists,
            "datasets": datasets,
            "layout": _layout(pconfig).to_plotly_json(),
            "pconfig": pconfig.__dict__,
            # Parameters to be toggled by switch buttons:
            "active_dataset_idx": 0,
            "p_active": pconfig.p_active,
            "l_active": pconfig.l_active,
        }

    return html


def _layout(pconfig) -> go.Layout:
    """
    Layout object for the bar plot.
    """
    layout: go.Layout = base_layout(pconfig)
    layout.update(
        {
            "barmode": pconfig.stacking,
            "hovermode": "y unified",
            "showlegend": True,
            "yaxis": dict(
                # the plot is "transposed", so yaxis corresponds to the horizontal axis
                title=dict(text=pconfig.xlab),
                showgrid=False,
                categoryorder="category descending",
                automargin=True,
            ),
            "xaxis": dict(
                title=dict(text=pconfig.ylab),
            ),
            "legend": dict(
                orientation="h",
                yanchor="top",
                y=-0.2,
                xanchor="center",
                x=0.5,
            ),
        }
    )
    return layout


def _datasets_to_interactive_imgs(
    datasets: List[List[Dict]],
    pconfig: PConfig,
    pid: str,
) -> str:
    html = '<div class="mqc_hcplot_plotgroup">'

    btn_tmpl = (
        f'<button class="{{cls}} btn btn-default btn-sm {{active}}" data-target="{pid}">{{label}}' f"</button> \n"
    )
    # Counts / Percentages / Log Switches
    if pconfig.add_pct_tab or pconfig.add_log_tab:
        # html += '<div class="btn-group hc_switch_group"> \n'
        if pconfig.add_pct_tab:
            html += btn_tmpl.format(
                active=pconfig.p_active,
                label=pconfig.p_label,
                cls="switch_percent",
            )
        if pconfig.add_log_tab:
            html += btn_tmpl.format(
                active=pconfig.l_active,
                label=pconfig.l_label,
                cls="switch_log10",
            )
        # html += "</div> "
        if len(datasets) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(datasets) > 1:
        html += '<div class="btn-group dataset_switch_group">\n'
        for k, ds in enumerate(datasets):
            active = "active" if k == 0 else ""
            dl: Dict[str, str] = pconfig.data_labels[k]
            name = dl["name"]
            ylab = f'data-ylab="{dl["ylab"]}"' if "ylab" in dl else ""
            ymax = f'data-ylab="{dl["ymax"]}"' if "ymax" in dl else ""
            html += f'<button class="btn btn-default btn-sm {active}" {ylab} {ymax} data-dataset_index="{k}" data-target="{pid}">{name}</button>\n'
        html += "</div>\n\n"

    # Plot HTML
    html += """<div class="hc-plot-wrapper"{height}>
        <div id="{id}" class="hc-plot not_rendered hc-bar-plot"></div>
    </div>""".format(
        id=pid,
        height=f' style="height:{pconfig.height}px"' if pconfig.height else "",
    )
    # Close wrapping div
    html += "</div>"
    return html


def _datasets_to_flat_imgs(
    datasets: List,
    samples_lists: List,
    pconfig: PConfig,
    pid: str,
    pids: List[str],
) -> str:
    """
    Create an HTML for a static plot version.
    """
    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )
    html += f'<div class="mqc_mplplot_plotgroup" id="{pid}">'

    # Counts / Percentages Switch
    if pconfig.add_pct_tab and not config.simple_output:
        html += f'<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_setcountspcnt"> \n\
            <button class="btn btn-default btn-sm {pconfig.c_active} counts">{pconfig.c_label}</button> \n\
            <button class="btn btn-default btn-sm {pconfig.p_active} pcnt">{pconfig.p_label}</button> \n\
        </div> '
        if len(datasets) > 1:
            html += " &nbsp; &nbsp; "

    # Buttons to cycle through different datasets
    if len(datasets) > 1 and not config.simple_output:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for pidx, ds in enumerate(datasets):
            pid = pids[pidx]
            active = "active" if pidx == 0 else ""
            name = pconfig.data_labels[pidx]["name"]
            html += f'<button class="btn btn-default btn-sm {active}" data-target="#{pid}">{name}</button>\n'
        html += "</div>\n\n"

        # Collect all categories, and fill in with zeroes for samples that having
        # any of cats missing
        cat_to_color = dict()
        for ds in datasets:
            for d in ds:
                cat_to_color[d["name"]] = d["color"]
        for data_by_cat in datasets:
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
        for data_by_cat in datasets:
            data_by_cat.sort(key=lambda x: x["name"])

    # Finally, build and save plots
    for pidx, (pid, samples, data_by_cat) in enumerate(zip(pids, samples_lists, datasets)):
        html += _dataset_to_imgs(pidx, pid, samples, data_by_cat, pconfig)

    # Close wrapping div
    html += "</div>"
    return html


def _save_data_file(pid, data_by_cat, samples):
    """
    Save plot data to file.
    """
    fdata = {}
    for d in data_by_cat:
        for didx, dval in enumerate(d["data"]):
            s_name = samples[didx]
            if s_name not in fdata:
                fdata[s_name] = dict()
            fdata[s_name][d["name"]] = dval
    util_functions.write_data_file(fdata, pid)


def _dataset_to_imgs(pidx, pid, samples, data_by_cat, pconfig) -> str:
    """
    Build a static images for different views of a dataset (counts, percentages, log scale),
    return an HTML wrapper.
    """
    # Save plot data to file
    if pconfig.save_data_file:
        _save_data_file(pid, data_by_cat, samples)

    # # Switch out NaN for 0s
    # for idx, d in enumerate(data_by_cat):
    #     data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

    # Calculate log10 values
    if pconfig.add_log_tab:
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

    views = [
        View(
            data_by_cat,
            active=not pconfig.p_active,
            suffix="",
            label=pconfig.c_label,
            xaxis_tickformat="",
        ),
    ]
    if pconfig.add_pct_tab:
        views.append(
            View(
                [
                    {
                        "data": d["data_pct"],
                        "name": d["name"],
                        "color": d["color"],
                    }
                    for d in data_by_cat
                ],
                active=pconfig.p_active,
                suffix="_pc",
                label=pconfig.p_label,
                xaxis_tickformat=".0%",
            )
        )
    if pconfig.add_log_tab:
        views.append(
            View(
                [
                    {
                        "data": d["data_log"],
                        "name": d["name"],
                        "color": d["color"],
                    }
                    for d in data_by_cat
                ],
                active=False,
                suffix="_log",
                label=pconfig.l_label,
                xaxis_tickformat="",
            )
        )

    html = ""
    for view in views:
        html += _view_to_img(
            view,
            pidx,
            f"{pid}{view.suffix}",
            samples,
            pconfig,
        )
    return html


def _view_to_img(
    view: View,
    pidx: int,
    pid: str,
    samples: List[str],
    pconfig: PConfig,
) -> str:
    """
    Build one static image, return an HTML wrapper.
    """
    pid = f"{pid}{view.suffix}"

    # Should this plot be hidden on report load?
    hide_div = ""
    if pidx > 0 or not view.active:
        hide_div = ' style="display:none;"'

    fig = go.Figure(layout=_layout(pconfig))
    for d in view.data:
        fig.add_trace(
            go.Bar(
                y=samples,
                x=d["data"],
                name=d["name"],
                orientation="h",
                marker=dict(color=d.get("color"), line=dict(width=0)),
            ),
        )
    if view.xaxis_tickformat:
        fig.update_layout({"xaxis": {"tickformat": view.xaxis_tickformat}})

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
        return f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="data:image/png;base64,{b64_img}" /></div>'

    # Link to the saved image
    else:
        plot_relpath = Path(config.plots_dir_name) / "png" / f"{pid}.png"
        plot_relpath.parent.mkdir(parents=True, exist_ok=True)
        fig.write_image(
            plot_relpath,
            format="png",
            width=fig.layout.width,
            height=fig.layout.height,
            scale=1,
        )
        return f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'
