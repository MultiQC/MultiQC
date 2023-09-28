"""Plotly bargraph functionality."""
import math
from typing import Dict, List, cast

import plotly.graph_objects as go

from multiqc.utils import mqc_colour


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

    fig = go.Figure()
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

    # NOTE: Only first plot for now, if multiple plots. To support multiple plots,
    # need to figure out a way to make the switch buttons reactive, so you can
    # switch both between the datasets and log/percentages/counts at the same time.
    data_by_cat = data_by_cat_lists[0]
    samples_lists = samples_lists[0]

    # Regular bar graph layout configurations here...
    fig.update_layout(
        cast(  # The function expects a dictionary, even though go.Layout works just fine
            dict,
            go.Layout(
                title=dict(
                    text=pconfig["title"],
                    xanchor="center",
                    x=0.5,
                    font=dict(size=20),
                ),
                yaxis=dict(
                    title=dict(text=pconfig.get("xlab")),
                    showgrid=False,
                    categoryorder="category descending",
                ),
                xaxis=dict(
                    title=dict(text=pconfig.get("ylab")),
                    gridcolor="rgba(0,0,0,0.1)",
                ),
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                hovermode="y unified",
                # template="simple_white",
                font={"color": "Black", "family": "Lucida Grande"},
                colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
            ),
        ),
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
        barmode="stack",
        height=plt_height,
    )

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
                    buttons=list(
                        [
                            dict(
                                args=["x", [data["data"] for data in data_by_cat]],
                                label=pconfig.get("cpswitch_counts_label", "Counts"),
                                method="restyle",
                            ),
                        ]
                        + (
                            [
                                dict(
                                    args=["x", pct_by_cat],
                                    label=pconfig.get("cpswitch_percent_label", "Percentages"),
                                    method="restyle",
                                ),
                            ]
                            if add_pct_tab
                            else []
                        )
                        + (
                            [
                                dict(
                                    args=["x", log_by_cat],
                                    label=pconfig.get("logswitch_label", "Log10"),
                                    method="restyle",
                                )
                            ]
                            if add_log_tab
                            else []
                        ),
                    ),
                    yanchor="top",
                    xanchor="left",
                    x=1,
                    y=1.2,
                ),
            ]
        )

    fig.add_annotation(
        xanchor="right",
        yanchor="bottom",
        x=1.1,
        y=-0.3,
        text="Created with MultiQC",
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=10, color="rgba(0,0,0,0.5)"),
    )

    img_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    return f'<div class="mqc_mplplot">{img_html}</div>'
