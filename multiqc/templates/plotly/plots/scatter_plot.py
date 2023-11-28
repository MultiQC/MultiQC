import base64
import io
import logging
import os
import random
import string
from pathlib import Path
from typing import Dict, List, Union, Optional

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots import get_template_mod
from multiqc.templates.plotly.plots.plot import View, PConfig, base_layout
from multiqc.utils import util_functions, config

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
SampleLineT = Dict[str, Union[str, List]]


class ScatterPlotConfig(PConfig):
    def __init__(self, pconfig: Dict, *args):
        super().__init__(pconfig, *args)
        self.tt_decimals: Optional[int] = pconfig.get("tt_decimals")
        self.tt_suffix: str = pconfig.get("tt_suffix", "")


def scatter_plot(
    datasets: List[List[SampleLineT]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    return add_to_report(
        datasets=datasets,
        pconfig=ScatterPlotConfig(pconfig, len(datasets)),
        report=report,
    )


def add_to_report(
    datasets: List[List[SampleLineT]],
    pconfig: ScatterPlotConfig,
    report,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
    is_static_suf = "static_" if config.plots_force_flat else ""
    pid = report.save_htmlid(pconfig.id)
    if pid is None:  # ID of the plot group
        pid = report.save_htmlid(f"mqc_{is_static_suf}plot_{uniq_suffix}")

    if config.plots_force_flat:
        return _datasets_to_flat_imgs(
            datasets,
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
            "plot_type": "xy_line",
            "datasets": datasets,
            "layout": _layout(pconfig).to_plotly_json(),
            "pconfig": pconfig.__dict__,
            # Parameters to be toggled by switch buttons:
            "active_dataset_idx": 0,
            "l_active": pconfig.l_active,
        }
        return html


def _layout(pconfig: ScatterPlotConfig) -> go.Layout:
    """
    Customise plot layout.
    """
    layout = base_layout(pconfig)

    layout.showlegend = False

    # Make a tooltip always show on hover over any point on plot
    layout.hoverdistance = -1
    # A tooltip will show numbers for all lines crossing this vertical line
    layout.hovermode = "x"
    # Default precision for floating numbers is too high - allowing to override it
    if pconfig.tt_decimals is not None:
        layout.yaxis.hoverformat = f".{pconfig.tt_decimals}f"

    return layout


def _datasets_to_interactive_imgs(
    datasets: List[List[SampleLineT]],
    pconfig: ScatterPlotConfig,
    pid: str,
) -> str:
    html = '<div class="mqc_hcplot_plotgroup">'
    btn_tmpl = (
        f'<button class="{{cls}} btn btn-default btn-sm {{active}}" data-target="{pid}">{{label}}' f"</button> \n"
    )
    # Counts / Percentages / Log Switches
    if pconfig.add_log_tab:
        # html += '<div class="btn-group hc_switch_group"> \n'
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
            xlab = f'data-ylab="{dl["xlab"]}"' if "xlab" in dl else ""
            html += (
                f'<button class="btn btn-default btn-sm {active}" {ylab} {ymax} {xlab} data-dataset_index="{k}" '
                f'data-target="{pid}">{name}</button>\n'
            )
        html += "</div>\n\n"

    # Plot HTML
    html += """<div class="hc-plot-wrapper"{height}>
        <div id="{id}" class="hc-plot not_rendered hc-bar-plot"></div>
    </div></div>""".format(
        id=pid,
        height=f' style="height:{pconfig.height}px"' if pconfig.height else "",
    )

    return html


def _datasets_to_flat_imgs(
    datasets: List[List[SampleLineT]],
    pconfig: ScatterPlotConfig,
    pid: str,
    pids: List[str],
    no_js: bool = False,
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

    # Buttons to cycle through different datasets
    if len(datasets) > 1 and not no_js:
        html += '<div class="btn-group mpl_switch_group mqc_mplplot_bargraph_switchds">\n'
        for k, p in enumerate(datasets):
            pid = pids[k]
            active = "active" if k == 0 else ""
            name = pconfig.data_labels[k]["name"]
            html += f'<button class="btn btn-default btn-sm {active}" data-target="#{pid}">{name}</button>\n'
        html += "</div>\n\n"

    # Go through datasets creating plots
    for pidx, (pid, dataset) in enumerate(zip(pids, datasets)):
        html += _dataset_to_imgs(pidx, pid, dataset, pconfig)

    # Close wrapping div
    html += "</div>"
    return html


def _save_data_file(
    pid: str,
    dataset: List[SampleLineT],
    pconfig: ScatterPlotConfig,
) -> None:
    fdata = dict()
    last_cats = None
    shared_cats = True
    for ds in dataset:
        fdata[ds["name"]] = dict()

        # Check to see if all categories are the same
        if len(ds["data"]) > 0 and isinstance(ds["data"][0], list):
            if last_cats is None:
                last_cats = [x[0] for x in ds["data"]]
            elif last_cats != [x[0] for x in ds["data"]]:
                shared_cats = False

        for i, x in enumerate(ds["data"]):
            if isinstance(x, list):
                fdata[ds["name"]][x[0]] = x[1]
            else:
                try:
                    fdata[ds["name"]][pconfig.categories[i]] = x
                except Exception:
                    fdata[ds["name"]][str(i)] = x

    # Custom tsv output if the x-axis varies
    if not shared_cats and config.data_format == "tsv":
        fout = ""
        for ds in dataset:
            fout += "\t" + "\t".join([str(x[0]) for x in ds["data"]])
            fout += "\n{}\t".format(ds["name"])
            fout += "\t".join([str(x[1]) for x in ds["data"]])
            fout += "\n"
        with io.open(os.path.join(config.data_dir, f"{pid}.txt"), "w", encoding="utf-8") as f:
            print(fout.encode("utf-8", "ignore").decode("utf-8"), file=f)
    else:
        util_functions.write_data_file(fdata, pid)


def _dataset_to_imgs(
    pidx: int,
    pid: str,
    dataset: List[SampleLineT],
    pconfig: ScatterPlotConfig,
) -> str:
    """
    Build a static images for different views of a dataset (counts, log scale),
    return an HTML wrapper.
    """
    # Save plot data to file
    if pconfig.save_data_file:
        _save_data_file(
            pid,
            dataset,
            pconfig,
        )

    html = ""

    # Calculate log10 values
    if pconfig.add_log_tab:
        for ds in dataset:
            ds["data_log"] = []
            for x, y in ds["data"]:
                if y == 0:
                    y = 0
                else:
                    y = math.log10(y)
                ds["data_log"].append([x, y])

    views = [
        View(
            dataset,
            active=not pconfig.l_active,
            suffix="",
            label=pconfig.c_label,
            xaxis_tickformat="",
        ),
    ]
    if pconfig.add_log_tab:
        views.append(
            View(
                [
                    {
                        "data": ds["data_log"],
                        "name": ds["name"],
                        "color": ds["color"],
                    }
                    for ds in dataset
                ],
                active=False,
                suffix="_log",
                label=pconfig.l_label,
                xaxis_tickformat="",
            )
        )

    for view in views:
        html += _view_to_img(
            view,
            pidx,
            f"{pid}{view.suffix}",
            pconfig,
            config,
        )

    return html


def _view_to_img(
    view: View,
    pidx: int,
    pid: str,
    pconfig: ScatterPlotConfig,
    config,
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
    for sdata in view.data:
        if len(sdata["data"]) > 0 and isinstance(sdata["data"][0], list):
            x = [x[0] for x in sdata["data"]]
            y = [x[1] for x in sdata["data"]]
        else:
            x = [x for x in range(len(sdata["data"]))]
            y = sdata["data"]

        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=sdata["name"],
                mode="lines+markers",
                marker=dict(size=5),
            )
        )

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
            width=900,
            height=fig.layout.height,
            scale=1,
        )
        return f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'
