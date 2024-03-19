import copy
import dataclasses
import logging
from typing import Dict, List, Union

import plotly.graph_objects as go

from multiqc.plots.plotly import determine_barplot_height
from multiqc.plots.plotly.plot import Plot, PlotType, BaseDataset
from multiqc.utils import util_functions

logger = logging.getLogger(__name__)


# Type of single box (matching one sample)
BoxT = List[Union[int, float]]


def plot(list_of_data_by_sample: List[Dict[str, BoxT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param list_of_data_by_sample: each dataset is a dict mapping samples to either:
        * list of data points,
        * list of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    :param pconfig: Plot configuration dictionary
    :return: HTML with JS, ready to be inserted into the page
    """
    p = BoxPlot(
        pconfig,
        list_of_data_by_sample,
        max_n_samples=max(len(d) for d in list_of_data_by_sample),
    )

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
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
        is_log=False,
        is_pct=False,
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


class BoxPlot(Plot):
    def __init__(
        self,
        pconfig: Dict,
        list_of_data_by_sample: List[Dict[str, BoxT]],
        max_n_samples: int,
    ):
        super().__init__(PlotType.BOX, pconfig, n_datasets=len(list_of_data_by_sample))

        self.datasets: List[Dataset] = [
            Dataset.create(ds, data_by_sample) for ds, data_by_sample in zip(self.datasets, list_of_data_by_sample)
        ]

        height = determine_barplot_height(max_n_samples)

        self.layout.update(
            height=height,
            showlegend=False,
            boxgroupgap=0.1,
            boxgap=0.2,
            colorway=[],  # no need to color code
            yaxis=dict(
                automargin=True,  # to make sure there is enough space for ticks labels
                categoryorder="trace",  # keep sample order
                hoverformat=self.layout.xaxis.hoverformat,
                ticksuffix=self.layout.xaxis.ticksuffix,
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                title=dict(text=self.layout.yaxis.title.text),
                hoverformat=self.layout.yaxis.hoverformat,
                ticksuffix=self.layout.yaxis.ticksuffix,
            ),
            hovermode="y",
            hoverlabel=dict(
                bgcolor="rgba(255, 255, 255, 0.8)",
                font=dict(color="black"),
            ),
        )

    def save_data_file(self, dataset: Dataset) -> None:
        vals_by_sample = {}
        for sample, values in zip(dataset.samples, dataset.data):
            vals_by_sample[sample] = values
        util_functions.write_data_file(vals_by_sample, dataset.uid)
