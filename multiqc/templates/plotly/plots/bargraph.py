"""Plotly bargraph functionality."""
import base64
import io
import math
import os
import string
from collections import namedtuple
from pathlib import Path
from random import random
from typing import Dict, List

import plotly.graph_objects as go

from multiqc.templates.plotly.plots import basic_figure
from multiqc.utils import config, report, util_functions

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


# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def _bargraph_plotly_plot(
    data_by_cat: List[Dict],
    sample_names: List[str],
    pconfig: Dict,
    plt_height: int = None,
):
    plt_height = plt_height or pconfig.get("height")
    fig = basic_figure(pconfig, plt_height)
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
    pconfig: Dict = None,
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
    if pconfig is None:
        pconfig = {}

    plt_height = pconfig.get("height")
    if not plt_height:
        max_n_samples = max(len(samples) for samples in samples_lists)
        # Height has a default, then adjusted by the number of samples
        plt_height = max_n_samples // 186  # Default, empirically determined
        plt_height = max(600, plt_height)  # At least 512px tall
        plt_height = min(2560, plt_height)  # Cap at 2560px tall

    uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
    is_static_suf = "static_" if config.plots_force_flat else ""
    # Plot group ID
    if pconfig.get("id") is None:
        pconfig["id"] = f"mqc_{is_static_suf}plot_{uniq_suffix}"
    # Sanitise plot ID and check for duplicates
    pconfig["id"] = report.save_htmlid(pconfig["id"])

    if config.plots_force_flat:
        return _static_bargraph(
            data_by_cat_lists,
            samples_lists,
            pconfig,
            plt_height,
        )

    fig = _bargraph_plotly_plot(
        data_by_cat_lists[0],
        samples_lists[0],
        pconfig,
        plt_height,
    )

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
            except:
                try:
                    name = pconfig["data_labels"][k]
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

    # for data, samples in zip(data_by_cat_lists, samples_lists):
    #     add_pct_tab = pconfig.get("cpswitch", True)
    #     add_log_tab = pconfig.get("logswitch", False)
    #     if add_pct_tab or add_log_tab:
    #         pct_by_cat = []
    #         if add_pct_tab:
    #             # Count totals for each category
    #             sums = [0 for _ in data[0]["data"]]
    #             for cat_idx, d in enumerate(data):
    #                 for sample_idx, v in enumerate(d["data"]):
    #                     if not math.isnan(v):
    #                         sums[sample_idx] += v
    #             # Now, calculate percentages for each category
    #             for cat_idx, d in enumerate(data):
    #                 values = [x for x in d["data"]]
    #                 if len(values) < len(samples):
    #                     values.extend([0] * (len(samples) - len(values)))
    #                 for key, var in enumerate(values):
    #                     sum_for_cat = sums[key]
    #                     if sum_for_cat == 0:
    #                         values[key] = 0
    #                     else:
    #                         values[key] = (float(var + 0.0) / float(sum_for_cat)) * 100
    #                 pct_by_cat.append(values)
    #
    #         log_by_cat = []
    #         if add_log_tab:
    #             for cat_idx, d in enumerate(data):
    #                 values = [x for x in d["data"]]
    #                 if len(values) < len(samples):
    #                     values.extend([0] * (len(samples) - len(values)))
    #                 for key, var in enumerate(values):
    #                     if var == 0:
    #                         values[key] = 0
    #                     else:
    #                         values[key] = math.log10(var)
    #                 log_by_cat.append(values)

    # fig.update_layout(
    #     updatemenus=[
    #         go.layout.Updatemenu(
    #             type="buttons",
    #             direction="left",
    #             buttons=[
    #                 b
    #                 for b in [
    #                     dict(
    #                         args=["x", [data["data"] for data in data]],
    #                         label=pconfig.get("cpswitch_counts_label", "Counts"),
    #                         method="restyle",
    #                     ),
    #                     dict(
    #                         args=["x", pct_by_cat],
    #                         label=pconfig.get("cpswitch_percent_label", "Percentages"),
    #                         method="restyle",
    #                     )
    #                     if add_pct_tab
    #                     else None,
    #                     dict(
    #                         args=["x", log_by_cat],
    #                         label=pconfig.get("logswitch_label", "Log10"),
    #                         method="restyle",
    #                     )
    #                     if add_log_tab
    #                     else None,
    #                 ]
    #                 if b
    #             ],
    #             yanchor="top",
    #             xanchor="left",
    #             x=1,
    #             y=1.2,
    #         ),
    #     ]
    # )

    # json = fig.to_plotly_json(full_html=False, include_plotlyjs=None)
    html = fig.to_html(full_html=False, include_plotlyjs=None)

    # html = """
    # <div class="plotly-plot-wrapper"{height}>
    #     <div id="{id}" class="plotly-plot not_rendered plotly-bar-plot"><small>loading..</small></div>
    # </div></div>""".format(
    #     id=pconfig["id"],
    #     height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    # )
    #
    # report.num_hc_plots += 1
    #
    # report.plot_data[pconfig["id"]] = {
    #     "plot_type": "bar_graph",
    #     "samples": samples_lists,
    #     "datasets": data_by_cat_lists,
    #     "config": pconfig,
    # }
    return html


def _static_bargraph(
    data_by_cat_lists: List,
    samples_lists: List,
    pconfig: Dict,
    plt_height: int = None,
):
    # Individual plot IDs
    pids = []
    for k in range(len(data_by_cat_lists)):
        try:
            name = pconfig["data_labels"][k]
        except KeyError:
            name = k + 1
        pid = report.save_htmlid(f"{pconfig['id']}_{name}", skiplint=True)
        pids.append(pid)

    html = (
        '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
        + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
        + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
    )
    html += f'<div class="mqc_mplplot_plotgroup" id="{pconfig["id"]}">'

    # Counts / Percentages Switch
    add_pct_tab = pconfig.get("cpswitch", True) is not False
    add_log_tab = pconfig.get("logswitch", False)
    pct_active = pconfig.get("cpswitch_c_active", True) is False or not add_pct_tab
    c_label = pconfig.get("cpswitch_counts_label", "Counts")
    p_label = pconfig.get("cpswitch_percent_label", "Percentages")

    if add_pct_tab and not config.simple_output:
        if not pct_active:
            c_active = "active"
            p_active = ""
        else:
            c_active = ""
            p_active = "active"
            pconfig["stacking"] = "percent"
        html += f'<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_setcountspcnt"> \n\
            <button class="btn btn-default btn-sm {c_active} counts">{c_label}</button> \n\
            <button class="btn btn-default btn-sm {p_active} pcnt">{p_label}</button> \n\
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
                name = pconfig["data_labels"][k]
            except:
                name = k + 1
            html += '<button class="btn btn-default btn-sm {a}" data-target="#{pid}">{n}</button>\n'.format(
                a=active, pid=pid, n=name
            )
        html += "</div>\n\n"

        # Collect all categories, and fill in with zeroes for samples that having any of cats missing
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
        if pconfig.get("save_data_file", True):
            util_functions.write_data_file(fdata, pids[pidx])

        # Switch out NaN for 0s
        for idx, d in enumerate(data_by_cat):
            data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        # Calculate percentages and log10 values
        pct_by_cat = []
        log_by_cat = []
        if add_pct_tab or add_log_tab:
            pct_by_cat = []
            if add_pct_tab:
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
                    pct_by_cat.append(values)

            log_by_cat = []
            if add_log_tab:
                for cat_idx, d in enumerate(data_by_cat):
                    values = [x for x in d["data"]]
                    if len(values) < len(samples):
                        values.extend([0] * (len(samples) - len(values)))
                    for key, var in enumerate(values):
                        if var == 0:
                            values[key] = 0
                        else:
                            values[key] = math.log10(var)
                    log_by_cat.append(values)

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
                active=not pct_active,
                suffix="",
                label=pconfig.get("cpswitch_counts_label", "Counts"),
            ),
        ]
        if add_pct_tab:
            views.append(
                View(
                    [
                        {
                            "data": vals,
                            "name": data_by_cat[i]["name"],
                            "color": data_by_cat[i]["color"],
                        }
                        for i, vals in enumerate(pct_by_cat)
                    ],
                    active=pct_active,
                    suffix="_pc",
                    label=pconfig.get("cpswitch_percent_label", "Percentages"),
                )
            )
        if add_log_tab:
            views.append(
                View(
                    [
                        {
                            "data": vals,
                            "name": data_by_cat[i]["name"],
                            "color": data_by_cat[i]["color"],
                        }
                        for i, vals in enumerate(log_by_cat)
                    ],
                    active=False,
                    suffix="_log",
                    label=pconfig.get("logswitch_label", "Log10"),
                )
            )

        for view in views:
            pid = f"{pids[pidx]}{view.suffix}"
            plot = _bargraph_plotly_plot(
                data_by_cat=view.values_by_cat,
                sample_names=samples_lists[pidx],
                pconfig=pconfig,
                plt_height=plt_height,
            )

            # Should this plot be hidden on report load?
            hide_div = ""
            if pidx > 0 or not view.active:
                hide_div = ' style="display:none;"'

            # Save the plot to the data directory if export is requested
            if config.export_plots:
                for fformat in config.export_plot_formats:
                    # Make the directory if it doesn't already exist
                    plot_dir = os.path.join(config.plots_dir, fformat)
                    if not os.path.exists(plot_dir):
                        os.makedirs(plot_dir)
                    # Save the plot
                    plot_fn = os.path.join(plot_dir, "{}.{}".format(pid, fformat))
                    plot.write_image(
                        plot_fn, format=fformat, width=plot.layout.width, height=plot.layout.height, scale=1
                    )

            # Output the figure to a base64 encoded string
            if getattr(get_template_mod(), "base64_plots", True) is True:
                img_buffer = io.BytesIO()
                plot.write_image(img_buffer, format="png", width=plot.layout.width, height=plot.layout.height, scale=1)
                b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
                img_buffer.close()
                html += '<div class="mqc_mplplot" id="{}"{}><img src="data:image/png;base64,{}" /></div>'.format(
                    pid, hide_div, b64_img
                )

            # Link to the saved image
            else:
                plot_relpath = Path(config.plots_dir_name) / "png" / f"{pid}.png"
                plot_relpath.parent.mkdir(parents=True, exist_ok=True)
                plot.write_image(
                    plot_relpath, format="png", width=plot.layout.width, height=plot.layout.height, scale=1
                )
                html += f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'

    # Close wrapping div
    html += "</div>"
    return html
