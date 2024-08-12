import copy
import logging
from typing import Dict, List, Union

import plotly.graph_objects as go  # type: ignore

from multiqc.plots.plotly import determine_barplot_height
from multiqc.plots.plotly.plot import PlotType, BaseDataset, Plot, PConfig
from multiqc import report

logger = logging.getLogger(__name__)


class BoxPlotConfig(PConfig):
    sort_samples: bool = True


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
        is_log=False,
        is_pct=False,
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
        vals_by_sample = {}
        for sample, values in zip(self.samples, self.data):
            vals_by_sample[sample] = values
        report.write_data_file(vals_by_sample, self.uid)


class BoxPlot(Plot[Dataset]):
    datasets: List[Dataset]

    @staticmethod
    def create(
        pconfig: BoxPlotConfig,
        list_of_data_by_sample: List[Dict[str, BoxT]],
    ) -> "BoxPlot":
        model = Plot.initialize(
            plot_type=PlotType.BOX,
            pconfig=pconfig,
            n_samples_per_dataset=[len(x) for x in list_of_data_by_sample],
        )

        model.datasets = [
            Dataset.create(ds, data_by_sample) for ds, data_by_sample in zip(model.datasets, list_of_data_by_sample)
        ]

        max_n_samples = max(len(x) for x in list_of_data_by_sample) if list_of_data_by_sample else 0
        height = determine_barplot_height(max_n_samples)

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
