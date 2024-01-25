"""Plotly bargraph functionality."""
import copy
import dataclasses
import logging
import re
from typing import Dict, List

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
from multiqc.utils import util_functions

logger = logging.getLogger(__name__)


def plot(
    cats_lists: List[List[Dict]],
    samples_lists: List[List[str]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param cats_lists: each dataset is a list of dicts with the keys: {name, color, data},
        where `name` is the category name, `color` is the color of the bar,
        and `data` is a list of values for each sample. Each outer list will
        correspond a separate tab.
    :param samples_lists: list of lists of bar names (that is, sample names). Similarly,
        each outer list will correspond to a separate tab.
    :param pconfig: Plot configuration dictionary
    :return: HTML with JS, ready to be inserted into the page
    """
    p = BarPlot(
        pconfig,
        cats_lists,
        samples_lists,
        max_n_samples=max([len(samples) for samples in samples_lists]),
    )

    from multiqc.utils import report

    return p.add_to_report(report)


def _split_long_string(s: str, max_width=80) -> List[str]:
    """
    Split string into lines of max_width characters
    """
    lines = []
    current_line = ""
    words = re.split(r"(\W+)", s)
    for word in words:
        if len(current_line + word) <= max_width:
            current_line += word
        else:
            if current_line:
                lines.append(current_line)
            current_line = word

    if current_line:
        lines.append(current_line)

    return lines


class BarPlot(Plot):
    @dataclasses.dataclass
    class Dataset(BaseDataset):
        cats: List[Dict]
        samples: List[str]

        @staticmethod
        def create(
            dataset: BaseDataset,
            cats: List[Dict],
            samples: List[str],
        ) -> "BarPlot.Dataset":
            for cat in cats:
                cat["name"] = "<br>".join(_split_long_string(cat["name"]))
            ds = BarPlot.Dataset(
                **dataset.__dict__,
                cats=cats,
                samples=samples,
            )
            return ds

    def __init__(self, pconfig: Dict, cats_lists: List, samples_lists: List, max_n_samples: int):
        # swap x and y axes parameters: the bar plot is "transposed", so yaxis corresponds to the horizontal axis
        for x_param in ["xmin", "xmax", "xCeiling"]:
            y_param = "y" + x_param[1:]
            pconfig[x_param], pconfig[y_param] = pconfig.get(y_param), pconfig.get(x_param)

        super().__init__(PlotType.BAR, pconfig, len(cats_lists))
        if len(cats_lists) != len(samples_lists):
            raise ValueError("Number of datasets and samples lists do not match")

        # Extend each dataset object with a list of samples
        self.datasets: List[BarPlot.Dataset] = [
            BarPlot.Dataset.create(d, cats=cats, samples=samples)
            for d, cats, samples in zip(self.datasets, cats_lists, samples_lists)
        ]

        # set height to be proportional to the number of samples
        MIN_PLOT_HEIGHT = 400
        MAX_PLOT_HEIGHT = 2560
        px_per_sample = 50
        if max_n_samples > 5:
            px_per_sample = 35
        if max_n_samples > 10:
            px_per_sample = 30
        if max_n_samples > 20:
            px_per_sample = 25
        if max_n_samples > 30:
            px_per_sample = 20
        height = max_n_samples * px_per_sample

        # set height to be proportional to the number of categories
        max_n_cats = max([len(dataset.cats) for dataset in self.datasets])
        height = min(800, max(height, 19 * max_n_cats + 140))

        height = min(MAX_PLOT_HEIGHT, height)
        height = max(MIN_PLOT_HEIGHT, height)
        self.layout.height = height

        # swap x and y axes: the bar plot is "transposed", so yaxis corresponds to the horizontal axis
        barmode = "stack" if self.p_active else "relative"
        if "stacking" in pconfig and pconfig["stacking"] == "percentage":
            barmode = "relative"
        self.layout.update(
            dict(
                showlegend=True,
                legend=None,  # default legend location
                barmode=barmode,
                hovermode="y unified",
                # xaxis_title_text=self.layout.xaxis.title.text,
                # yaxis_title_text=self.layout.yaxis.title.text,
                yaxis=dict(
                    showgrid=False,
                    categoryorder="category descending",  # otherwise the bars will be in reversed order to sample order
                    automargin=True,  # to make sure there is enough space for ticks labels
                    title=None,
                ),
                xaxis=dict(
                    title=dict(text=self.layout.yaxis.title.text),
                    range=[0, self.layout.xaxis.range[1]],
                ),
            )
        )

        self.trace_params = dict(
            orientation="h",
            marker=dict(line=dict(width=0)),
        )

        # Expand data with zeroes if there are fewer values than samples
        for dataset in self.datasets:
            for cat in dataset.cats:
                if len(cat["data"]) < len(dataset.samples):
                    cat["data"].extend([0] * (len(dataset.samples) - len(cat["data"])))

        # Calculate and save percentages
        if self.add_pct_tab:
            for pidx, dataset in enumerate(self.datasets):
                # Count totals for each category
                sums = [0 for _ in dataset.cats[0]["data"]]
                for cat in dataset.cats:
                    for sample_idx, val in enumerate(cat["data"]):
                        if not math.isnan(val):
                            sums[sample_idx] += val

                # Now, calculate percentages for each category
                for cat in dataset.cats:
                    values = [x for x in cat["data"]]
                    for key, var in enumerate(values):
                        sum_for_cat = sums[key]
                        if sum_for_cat == 0:
                            values[key] = 0
                        else:
                            values[key] = float(var + 0.0) / float(sum_for_cat)
                    cat["data_pct"] = values

        if self.add_log_tab:
            # Sorting from small to large so the log switch makes sense
            for dataset in self.datasets:
                dataset.cats.sort(key=lambda x: sum(x["data"]))

    @staticmethod
    def axis_controlled_by_switches() -> List[str]:
        """
        Return a list of axis names that are controlled by the log10 scale and percentage
        switch buttons
        """
        return ["xaxis"]

    def create_figure(
        self,
        layout: go.Layout,
        dataset: Dataset,
        is_log=False,
        is_pct=False,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)
        for cat in dataset.cats:
            data = cat["data"]
            if is_pct:
                data = cat["data_pct"]

            params = copy.deepcopy(self.trace_params)
            marker = params["marker"]
            color = cat.get("color")
            if color:
                marker["color"] = color

            fig.add_trace(
                go.Bar(
                    y=dataset.samples,
                    x=data,
                    name=cat["name"],
                    **params,
                ),
            )
        return fig

    def save_data_file(self, dataset: Dataset) -> None:
        fdata = {}
        for cat in dataset.cats:
            for d_idx, d_val in enumerate(cat["data"]):
                s_name = dataset.samples[d_idx]
                if s_name not in fdata:
                    fdata[s_name] = dict()
                fdata[s_name][cat["name"]] = d_val
        util_functions.write_data_file(fdata, dataset.uid)
