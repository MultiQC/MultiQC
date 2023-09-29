"""Plotly bargraph functionality."""
import math
from typing import Dict, List

import plotly.graph_objects as go

from multiqc.templates.plotly.plots import basic_figure
from multiqc.utils import report


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
    max_n_samples = max(len(samples) for samples in samples_lists)
    # Height has a default, then adjusted by the number of samples
    plt_height = max_n_samples // 186  # Default, empirically determined
    plt_height = max(600, plt_height)  # At least 512px tall
    plt_height = min(2560, plt_height)  # Cap at 2560 tall

    fig = basic_figure(pconfig, plt_height)

    for data in data_by_cat_lists[0]:
        fig.add_trace(
            go.Bar(
                y=samples_lists[0],
                x=data["data"],
                name=data["name"],
                orientation="h",
                marker=dict(color=data.get("color"), line=dict(width=0)),
            ),
        )

    for data_by_cat, samples in zip(data_by_cat_lists, samples_lists):
        add_pct_tab = pconfig.get("cpswitch", True)
        add_log_tab = pconfig.get("logswitch", False)
        if add_pct_tab or add_log_tab:
            pct_by_cat = []
            if add_pct_tab:
                # Count totals for each category
                sums = [0 for _ in data_by_cat[0]["data"]]
                for cat_idx, d in enumerate(data_by_cat):
                    for sample_idx, v in enumerate(d["data"]):
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

            fig.update_layout(
                updatemenus=[
                    go.layout.Updatemenu(
                        type="buttons",
                        direction="left",
                        buttons=[
                            dict(
                                args=["x", [data["data"] for data in data_by_cat]],
                                label=pconfig.get("cpswitch_counts_label", "Counts"),
                                method="restyle",
                            ),
                            dict(
                                args=["x", pct_by_cat],
                                label=pconfig.get("cpswitch_percent_label", "Percentages"),
                                method="restyle",
                            )
                            if add_pct_tab
                            else None,
                            dict(
                                args=["x", log_by_cat],
                                label=pconfig.get("logswitch_label", "Log10"),
                                method="restyle",
                            )
                            if add_log_tab
                            else None,
                        ],
                        yanchor="top",
                        xanchor="left",
                        x=1,
                        y=1.2,
                    ),
                ]
            )

    # json = fig.to_plotly_json(full_html=False, include_plotlyjs=None)
    html = fig.to_html(full_html=False, include_plotlyjs=None)
    html = """
    <div class="plotly-plot-wrapper"{height}>
        <div id="{id}" class="plotly-plot not_rendered plotly-bar-plot"><small>loading..</small></div>
    </div></div>""".format(
        id=pconfig["id"],
        height=f' style="height:{pconfig["height"]}px"' if "height" in pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[pconfig["id"]] = {
        "plot_type": "bar_graph",
        "samples": samples_lists,
        "datasets": data_by_cat_lists,
        "config": pconfig,
    }
    return html
