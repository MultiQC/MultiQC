"""Plotly bargraph functionality."""
import logging
from typing import Dict, List

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot
from multiqc.utils import config, util_functions

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

    # Sorting from small to large so the log switch makes sense
    for categories in datasets:
        categories.sort(key=lambda x: sum(x["data"]))

    for samples, categories in zip(samples_lists, datasets):
        for cat in categories:
            cat["samples"] = samples
            # Extend with zeroes if there are fewer values than samples
            if len(cat["data"]) < len(cat["samples"]):
                cat["data"].extend([0] * (len(cat["samples"]) - len(cat["data"])))

    p = BarPlot(pconfig, len(datasets), max(len(samples) for samples in samples_lists))

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

    return p.add_to_report(datasets, report)


class BarPlot(Plot):
    def __init__(self, pconfig: Dict, max_n_samples: int, *args):
        super().__init__("bar_graph", pconfig, *args)

        self.add_pct_tab = pconfig.get("cpswitch", True) is not False
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.p_active = ""
        self.stacking = pconfig.get("stacking", "stack")
        if self.add_pct_tab and pconfig.get("cpswitch_c_active", True) is not True and not config.simple_output:
            self.c_active = ""
            self.p_active = "active"
            self.stacking = "relative"

        if not self.height:
            # Height has a default, then adjusted by the number of samples
            self.height = max_n_samples // 186  # Default, empirically determined
            self.height = max(600, self.height)
            self.height = min(2560, self.height)

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()
        layout.update(
            {
                "barmode": self.stacking,
                "hovermode": "y unified",
                "showlegend": True,
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

    # def flat_plot(
    #     self,
    #     datasets: List,
    #     pid: str,
    #     pids: List[str],
    # ) -> str:
    #     """
    #     Create an HTML for a static plot version.
    #     """
    #     html = super().flat_plot(datasets, pid, pids)
    #
    #     # Buttons to cycle through different datasets
    #     if len(datasets) > 1 and not config.simple_output:
    #         # Collect all categories, and fill in with zeroes for samples that having
    #         # any of cats missing
    #         cat_to_color = dict()
    #         for ds in datasets:
    #             for d in ds:
    #                 cat_to_color[d["name"]] = d["color"]
    #         for data_by_cat in datasets:
    #             for cat, color in cat_to_color.items():
    #                 samples = data_by_cat[0]["samples"]
    #                 if cat not in [d["name"] for d in data_by_cat]:
    #                     data_by_cat.append(
    #                         {
    #                             "name": cat,
    #                             "color": color,
    #                             "data": [0 for _ in samples],
    #                             "data_pct": [0.0 for _ in samples],
    #                             "data_log": [0.0 for _ in samples],
    #                         }
    #                     )
    #         # Sort categories by name
    #         for data_by_cat in datasets:
    #             data_by_cat.sort(key=lambda x: x["name"])
    #
    #     # Finally, build and save plots
    #     for pidx, (pid, data_by_cat) in enumerate(zip(pids, datasets)):
    #         html += self.dataset_to_imgs(pidx, pid, data_by_cat)
    #
    #     # Close wrapping div
    #     html += "</div>"
    #     return html

    def populate_figure(self, fig: go.Figure, dataset: List[Dict], is_log=False, is_pct=False) -> go.Figure:
        fig = super()._make_fig(dataset, is_log, is_pct)

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
                    marker=dict(color=cat.get("color"), line=dict(width=0)),
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
