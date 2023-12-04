"""Plotly bargraph functionality."""
import logging
from typing import Dict, List

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType
from multiqc.utils import util_functions

logger = logging.getLogger(__name__)


def plot(
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
    :return: HTML with JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    for samples, categories in zip(samples_lists, datasets):
        for cat in categories:
            cat["samples"] = samples
            # Extend with zeroes if there are fewer values than samples
            if len(cat["data"]) < len(cat["samples"]):
                cat["data"].extend([0] * (len(cat["samples"]) - len(cat["data"])))

    p = BarPlot(
        pconfig,
        num_datasets=len(datasets),
        max_n_samples=max([len(samples) for samples in samples_lists]),
    )

    # Calculate and save percentages
    if p.add_pct_tab:
        for pidx, categories in enumerate(datasets):
            # Count totals for each category
            sums = [0 for _ in categories[0]["data"]]
            for cat in categories:
                for sample_idx, val in enumerate(cat["data"]):
                    if not math.isnan(val):
                        sums[sample_idx] += val

            # Now, calculate percentages for each category
            for cat in categories:
                values = [x for x in cat["data"]]
                for key, var in enumerate(values):
                    sum_for_cat = sums[key]
                    if sum_for_cat == 0:
                        values[key] = 0
                    else:
                        values[key] = float(var + 0.0) / float(sum_for_cat)
                cat["data_pct"] = values

    if p.add_log_tab:
        # Sorting from small to large so the log switch makes sense
        for categories in datasets:
            categories.sort(key=lambda x: sum(x["data"]))

    return p.add_to_report(datasets, report)


class BarPlot(Plot):
    def __init__(self, pconfig: Dict, max_n_samples: int, num_datasets: int):
        super().__init__(PlotType.BAR, pconfig, num_datasets=num_datasets)

        if not self.height:
            # Height has a default, then adjusted by the number of samples
            self.height = max_n_samples * 10
            self.height = max(500, self.height)
            self.height = min(2560, self.height)

        self.stacking = pconfig.get("stacking", "stack" if self.p_active else "relative")

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()
        layout.update(
            {
                "barmode": self.stacking,
                "hovermode": "y unified",
                "yaxis": dict(
                    # the plot is "transposed", so yaxis corresponds to the horizontal axis
                    title=dict(text=self.xlab),
                    showgrid=False,
                    categoryorder="category descending",
                    automargin=True,
                ),
                "xaxis": dict(
                    title=dict(text=self.ylab),
                ),
            }
        )
        return layout

    def populate_figure(self, fig: go.Figure, dataset: List[Dict], is_log=False, is_pct=False) -> go.Figure:
        categories: List[Dict] = dataset
        for cat in categories:
            data = cat["data"]
            if is_pct:
                data = cat["data_pct"]

            fig.add_trace(
                go.Bar(
                    y=cat["samples"],
                    x=data,
                    name=cat["name"],
                    orientation="h",
                    marker=dict(
                        color=cat.get("color"),
                        line=dict(width=0),
                    ),
                ),
            )
        return fig

    def save_data_file(self, dataset: List[Dict], uid: str) -> None:
        fdata = {}
        for d in dataset:
            for d_idx, d_val in enumerate(d["data"]):
                s_name = d["samples"][d_idx]
                if s_name not in fdata:
                    fdata[s_name] = dict()
                fdata[s_name][d["name"]] = d_val
        util_functions.write_data_file(fdata, uid)
