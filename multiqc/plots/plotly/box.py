import copy
import logging
from typing import Dict, List, Optional, Tuple, Union

import plotly.graph_objects as go  # type: ignore

from multiqc import report
from multiqc.plots.plotly import determine_barplot_height
from multiqc.plots.plotly.plot import BaseDataset, PConfig, Plot, PlotType
from multiqc.types import SampleName

logger = logging.getLogger(__name__)


class BoxPlotConfig(PConfig):
    sort_samples: bool = True

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("boxplot",), **data)


# Type of single box (matching one sample)
BoxT = List[Union[int, float]]


def plot(
    list_of_data_by_sample: List[Dict[str, BoxT]],
    pconfig: BoxPlotConfig,
) -> "BoxPlot":
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param list_of_data_by_sample: each dataset is a dict mapping samples to either:
        * list of data points,
        * list of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    :param pconfig: Plot configuration dictionary
    :return: HTML with JS, ready to be inserted into the page
    """
    return BoxPlot.create(
        pconfig=pconfig,
        list_of_data_by_sample=list_of_data_by_sample,
    )


class Dataset(BaseDataset):
    data: List[BoxT]
    samples: List[str]

    @staticmethod
    def create(
        dataset: BaseDataset,
        data_by_sample: Dict[str, BoxT],
    ) -> "Dataset":
        dataset = Dataset(
            **dataset.__dict__,
            data=list(data_by_sample.values()),
            samples=list(data_by_sample.keys()),
        )
        # Need to reverse samples as the box plot will show them reversed
        dataset.samples = list(reversed(dataset.samples))
        dataset.data = list(reversed(dataset.data))

        dataset.trace_params.update(
            boxpoints="outliers",
            jitter=0.5,
            orientation="h",
            marker=dict(
                color="#4899e8",  # just use blue to indicate interactivity
            ),
            # to remove the redundant sample name before "median" in the unified hover box
            hoverinfo="x",
        )
        dataset.layout["yaxis"]["title"] = None
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)

        for sname, values in zip(self.samples, self.data):
            params = copy.deepcopy(self.trace_params)
            fig.add_trace(
                go.Box(
                    x=values,
                    name=sname,
                    **params,
                ),
            )
        return fig

    def save_data_file(self) -> None:
        vals_by_sample: Dict[str, BoxT] = {}
        for sample, values in zip(self.samples, self.data):
            vals_by_sample[sample] = values
        report.write_data_file(vals_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        """Format dataset as a markdown table with basic statistics"""
        prompt = "|Sample|Min|Q1|Median|Q3|Max|Mean|\n"
        prompt += "|---|---|---|---|---|---|---|\n"

        suffix = ""
        if self.layout["xaxis"]["ticksuffix"]:
            suffix = " " + self.layout["xaxis"]["ticksuffix"]

        for sample, values in zip(self.samples, self.data):
            # Skip samples with no data
            if len(values) == 0:
                continue

            # Use pseudonym if available, otherwise use original sample name
            pseudonym = report.anonymize_sample_name(sample)

            # Calculate statistics
            sorted_vals = sorted(values)
            n = len(sorted_vals)

            min_val = sorted_vals[0]
            max_val = sorted_vals[-1]
            median = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
            q1 = sorted_vals[n // 4] if n >= 4 else sorted_vals[0]
            q3 = sorted_vals[3 * n // 4] if n >= 4 else sorted_vals[-1]
            mean = sum(values) / len(values)

            prompt += (
                f"|{pseudonym}|"
                f"{self.fmt_value_for_llm(min_val)}{suffix}|"
                f"{self.fmt_value_for_llm(q1)}{suffix}|"
                f"{self.fmt_value_for_llm(median)}{suffix}|"
                f"{self.fmt_value_for_llm(q3)}{suffix}|"
                f"{self.fmt_value_for_llm(max_val)}{suffix}|"
                f"{self.fmt_value_for_llm(mean)}{suffix}|\n"
            )
        return prompt


class BoxPlot(Plot[Dataset, BoxPlotConfig]):
    datasets: List[Dataset]

    def samples_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        for ds in self.datasets:
            names.extend(SampleName(sample) for sample in ds.samples)
        return names

    @staticmethod
    def create(
        pconfig: BoxPlotConfig,
        list_of_data_by_sample: List[Dict[str, BoxT]],
    ) -> "BoxPlot":
        model: Plot[Dataset, BoxPlotConfig] = Plot.initialize(
            plot_type=PlotType.BOX,
            pconfig=pconfig,
            n_samples_per_dataset=[len(x) for x in list_of_data_by_sample],
        )

        model.datasets = [
            Dataset.create(ds, data_by_sample) for ds, data_by_sample in zip(model.datasets, list_of_data_by_sample)
        ]

        max_n_samples = max(len(x) for x in list_of_data_by_sample) if list_of_data_by_sample else 0
        height: int = determine_barplot_height(max_n_samples)

        model.layout.update(
            height=height,
            showlegend=False,
            boxgroupgap=0.1,
            boxgap=0.2,
            colorway=[],  # no need to color code
            yaxis=dict(
                automargin=True,  # to make sure there is enough space for ticks labels
                categoryorder="trace",  # keep sample order
                hoverformat=model.layout.xaxis.hoverformat,
                ticksuffix=model.layout.xaxis.ticksuffix,
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                title=dict(text=model.layout.yaxis.title.text),
                hoverformat=model.layout.yaxis.hoverformat,
                ticksuffix=model.layout.yaxis.ticksuffix,
            ),
            hovermode="y",
            hoverlabel=dict(
                bgcolor="rgba(255, 255, 255, 0.8)",
                font=dict(color="black"),
            ),
        )
        return BoxPlot(**model.__dict__)
