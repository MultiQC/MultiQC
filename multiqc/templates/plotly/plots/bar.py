"""Plotly bargraph functionality."""
import copy
import dataclasses
import logging
import re
from collections import defaultdict
from typing import Dict, List

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
                cat["color"] = ",".join([f"{x:.2f}" for x in color.rgb])

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
            layout: go.Layout,
            is_log=False,
            is_pct=False,
        ) -> go.Figure:
            """
            Create a Plotly figure for a dataset
            """
            fig = go.Figure(layout=layout)
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

        # Set height to be proportional to the number of samples
        PADDING = 140  # plot header and footer
        height_per_bar = 35
        if max_n_samples > 5:
            height_per_bar = 30
        if max_n_samples > 10:
            height_per_bar = 25
        if max_n_samples > 20:
            height_per_bar = 20
        if max_n_samples > 30:
            height_per_bar = 15
        height = max_n_samples * height_per_bar + PADDING

        # Set height to also be proportional to the number of cats to fit a legend
        HEIGHT_PER_LEGEND_ITEM = 19
        n_cats = max([len(dataset.cats) for dataset in self.datasets])
        legend_height = HEIGHT_PER_LEGEND_ITEM * n_cats + PADDING
        # expand the plot to fit the legend:
        height = max(height, max(height, legend_height))
        # but not too much - if there are only 2 samples, we don't want the plot to be too high:
        height = min(800, max(height, legend_height))

        # now, limit the max and min height (plotly will start to automatically skip
        # some of the ticks on the left when the plot is too high)
        MIN_HEIGHT = 300
        MAX_HEIGHT = 2560
        height = min(MAX_HEIGHT, height)
        height = max(MIN_HEIGHT, height)

        # Set the barmode
        barmode = "relative"  # stacking, but drawing negative values below zero
        if "stacking" in pconfig and (pconfig["stacking"] in ["group", "normal", None]):
            barmode = "group"  # side by side

        self.layout.update(
            height=height,
            showlegend=True,
            barmode=barmode,
            yaxis=dict(
                showgrid=False,
                categoryorder="category descending",  # otherwise the bars will be in reversed order to sample order
                automargin=True,  # to make sure there is enough space for ticks labels
                title=None,
                hoverformat=self.layout.xaxis.hoverformat,
                ticksuffix=self.layout.xaxis.ticksuffix,
            ),
            xaxis=dict(
                title=dict(text=self.layout.yaxis.title.text),
                hoverformat=self.layout.yaxis.hoverformat,
                ticksuffix=self.layout.yaxis.ticksuffix,
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

        for dataset in self.datasets:
            if barmode == "group":
                # max category
                xmax_cnt = max(max(cat["data"][i] for cat in dataset.cats) for i in range(len(dataset.samples)))
                xmin_cnt = min(min(cat["data"][i] for cat in dataset.cats) for i in range(len(dataset.samples)))
            else:
                # max sum of all categories across all samples
                xmax_cnt = max(
                    sum(cat["data"][i] if cat["data"][i] > 0 else 0 for cat in dataset.cats)
                    for i in range(len(dataset.samples))
                )
                xmin_cnt = min(
                    sum(cat["data"][i] if cat["data"][i] < 0 else 0 for cat in dataset.cats)
                    for i in range(len(dataset.samples))
                )

            xmin_cnt = self.pconfig.get("ymin", xmin_cnt)
            xmax_cnt = self.pconfig.get("ymax", xmax_cnt)
            dataset.layout.update(
                yaxis=dict(
                    title=None,
                    hoverformat=dataset.layout["xaxis"]["hoverformat"],
                    ticksuffix=dataset.layout["xaxis"]["ticksuffix"],
                    # Prevent JavaScript from automatically parsing categorical values as numbers:
                    type="category",
                ),
                xaxis=dict(
                    title=dict(text=dataset.layout["yaxis"]["title"]["text"]),
                    hoverformat=dataset.layout["yaxis"]["hoverformat"],
                    ticksuffix=dataset.layout["yaxis"]["ticksuffix"],
                    range=[xmin_cnt, xmax_cnt],
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

            if dataset.layout["xaxis"]["hoverformat"] is None:
                if all(all(isinstance(x, float) or math.isnan(x) for x in cat["data"]) for cat in dataset.cats):
                    dataset.layout["xaxis"]["hoverformat"] = ",.2f"
                elif all(all(isinstance(x, int) or math.isnan(x) for x in cat["data"]) for cat in dataset.cats):
                    dataset.layout["xaxis"]["hoverformat"] = ",.0f"

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
                            sums[sample_idx] += abs(val)

                # Now, calculate percentages for each category
                for cat in dataset.cats:
                    values = [x for x in cat["data"]]
                    for sample_idx, val in enumerate(values):
                        sum_for_sample = sums[sample_idx]
                        if sum_for_sample == 0:
                            values[sample_idx] = 0
                        else:
                            values[sample_idx] = float(val + 0.0) / float(sum_for_sample) * 100.0
                    cat["data_pct"] = values

                if barmode == "group":
                    # calculating the min percentage range as well because it will be negative for negative values
                    xmin_pct = min(min(cat["data_pct"][i] for cat in dataset.cats) for i in range(len(dataset.samples)))
                else:
                    xmin_pct = min(
                        sum(cat["data_pct"][i] if cat["data_pct"][i] < 0 else 0 for cat in dataset.cats)
                        for i in range(len(dataset.samples))
                    )
                dataset.pct_range = [xmin_pct, 100]

        if self.add_log_tab:
            # Sorting from small to large so the log switch makes sense
            for dataset in self.datasets:
                dataset.cats.sort(key=lambda x: sum(x["data"]))
                # But reversing the legend so the largest bars are still on the top
                self.layout.legend.traceorder = "reversed"

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
