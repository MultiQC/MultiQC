"""Plotly bargraph functionality."""
import copy
import dataclasses
import logging
import re
from collections import defaultdict
from typing import Dict, List, Optional

import math
import plotly.graph_objects as go
import spectra

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
            # Post-process categories
            for cat in cats:
                # Split long category names
                if "name" not in cat:
                    raise ValueError(f"Bar plot {dataset.plot.id}: missing 'name' key in category")
                cat["name"] = "<br>".join(_split_long_string(cat["name"]))

                # Reformat color to be ready to add alpha in Plotly-JS
                color = spectra.html(cat["color"])
                cat["color"] = ",".join([str(x) for x in color.rgb])

                # Check that the number of samples is the same for all categories
                assert len(samples) == len(cat["data"])

            dataset = BarPlot.Dataset(
                **dataset.__dict__,
                cats=cats,
                samples=samples,
            )

            return dataset

        def create_figure(
            self,
            layout: Optional[go.Layout] = None,
            is_log=False,
            is_pct=False,
        ) -> go.Figure:
            """
            Create a Plotly figure for a dataset
            """
            fig = go.Figure(layout=layout or self.layout)
            for cat in self.cats:
                data = cat["data"]
                if is_pct:
                    data = cat["data_pct"]

                params = copy.deepcopy(self.trace_params)
                marker = params["marker"]
                marker["color"] = f"rgb({cat['color']})"

                fig.add_trace(
                    go.Bar(
                        y=self.samples,
                        x=data,
                        name=cat["name"],
                        **params,
                    ),
                )

            return fig

    def __init__(self, pconfig: Dict, cats_lists: List, samples_lists: List, max_n_samples: int):
        super().__init__(PlotType.BAR, pconfig, len(cats_lists))
        if len(cats_lists) != len(samples_lists):
            raise ValueError("Number of datasets and samples lists do not match")

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

        barmode = "stack"
        if "stacking" in pconfig and pconfig["stacking"] != "stack":
            barmode = "group"

        for dataset in self.datasets:
            xmax = self.pconfig.get("ymax")
            if xmax is None:
                if barmode == "stack":
                    xmax = max(sum(cat["data"][i] for cat in dataset.cats) for i in range(len(dataset.samples)))
                else:
                    xmax = max(max(cat["data"][i] for cat in dataset.cats) for i in range(len(dataset.samples)))
            dataset.layout.update(
                height=height,
                showlegend=True,
                barmode=barmode,
                yaxis=dict(
                    showgrid=False,
                    categoryorder="category descending",  # otherwise the bars will be in reversed order to sample order
                    automargin=True,  # to make sure there is enough space for ticks labels
                    title=None,
                    hoverformat=dataset.layout.xaxis.hoverformat,
                    ticksuffix=dataset.layout.xaxis.ticksuffix,
                ),
                xaxis=dict(
                    title=dict(text=dataset.layout.yaxis.title.text),
                    hoverformat=dataset.layout.yaxis.hoverformat,
                    ticksuffix=dataset.layout.yaxis.ticksuffix,
                    range=[0, xmax],
                ),
                # Re-initiate legend to reset to default legend location on the top right
                legend=go.layout.Legend(
                    # We use legend groups with subplots to simulate standard legend interactivity
                    # like we had a standard bar graph without subplots. We need to remove the space
                    # between the legend groups to make it look like a single legend.
                    tracegroupgap=0,
                ),
                hovermode="y unified",
                hoverlabel=dict(
                    bgcolor="rgba(255, 255, 255, 0.8)",
                    font=dict(color="black"),
                ),
            )
            dataset.trace_params.update(
                orientation="h",
                marker=dict(line=dict(width=0)),
                textposition="inside",
                insidetextanchor="start",
            )
            if "hovertemplate" in dataset.trace_params:
                # %{text} doesn't work for unified hovermode:
                dataset.trace_params["hovertemplate"] = dataset.trace_params["hovertemplate"].replace("%{text}", "")

            if dataset.layout.xaxis.hoverformat is None:
                if all(all(isinstance(x, float) for x in cat["data"]) for cat in dataset.cats):
                    dataset.layout.xaxis.hoverformat = ",.2f"
                elif all(all(isinstance(x, int) for x in cat["data"]) for cat in dataset.cats):
                    dataset.layout.xaxis.hoverformat = ",.0f"

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
                            values[key] = float(var + 0.0) / float(sum_for_cat) * 100.0
                    cat["data_pct"] = values

        if self.add_log_tab:
            # Sorting from small to large so the log switch makes sense
            for dataset in self.datasets:
                dataset.cats.sort(key=lambda x: sum(x["data"]))
                # But reversing the legend so the largest bars are still on the top
                dataset.layout.legend.traceorder = "reversed"

    @staticmethod
    def axis_controlled_by_switches() -> List[str]:
        """
        Return a list of axis names that are controlled by the log10 scale and percentage
        switch buttons
        """
        return ["xaxis"]

    def save_data_file(self, dataset: Dataset) -> None:
        val_by_cat_by_sample = defaultdict(dict)
        for cat in dataset.cats:
            for d_idx, d_val in enumerate(cat["data"]):
                s_name = dataset.samples[d_idx]
                val_by_cat_by_sample[s_name][cat["name"]] = d_val
        util_functions.write_data_file(val_by_cat_by_sample, dataset.uid)

    @staticmethod
    def tt_label() -> str:
        """Default tooltip label"""
        return "%{meta}: <b>%{x}</b>"
