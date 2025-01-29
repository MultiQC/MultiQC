"""MultiQC functions to plot a box plot"""

import copy
import logging
from typing import Any, Dict, List, Optional, OrderedDict, Tuple, Union, cast

import plotly.graph_objects as go  # type: ignore

from multiqc import report
from multiqc.plots.plot import BaseDataset, PConfig, Plot, PlotType, plot_anchor
from multiqc.plots.utils import determine_barplot_height
from multiqc.types import Anchor

logger = logging.getLogger(__name__)


class BoxPlotConfig(PConfig):
    sort_samples: bool = True

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("boxplot",), **data)


# Type of single box (matching one sample)
BoxT = List[Union[int, float]]


def plot(
    list_of_data_by_sample: Union[Dict[str, BoxT], List[Dict[str, BoxT]]],
    pconfig: Union[Dict[str, Any], BoxPlotConfig, None],
) -> "BoxPlot":
    """
    Plot a box plot. Expects either:
    - a dict mapping sample names to data point lists or dicts,
    - a dict mapping sample names to a dict of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    """
    pconf: BoxPlotConfig = cast(BoxPlotConfig, BoxPlotConfig.from_pconfig_dict(pconfig))

    anchor = plot_anchor(pconf)

    # Given one dataset - turn it into a list
    if not isinstance(list_of_data_by_sample, list):
        list_of_data_by_sample = [list_of_data_by_sample]

    for i in range(len(list_of_data_by_sample)):
        if isinstance(list_of_data_by_sample[0], OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconf.sort_samples:
            samples = sorted(list(list_of_data_by_sample[0].keys()))
            list_of_data_by_sample[i] = {s: list_of_data_by_sample[i][s] for s in samples}

    return BoxPlot.create(
        list_of_data_by_sample=list_of_data_by_sample,
        pconfig=pconf,
        anchor=anchor,
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


class BoxPlot(Plot[Dataset, BoxPlotConfig]):
    datasets: List[Dataset]

    @staticmethod
    def create(
        list_of_data_by_sample: List[Dict[str, BoxT]],
        pconfig: BoxPlotConfig,
        anchor: Anchor,
    ) -> "BoxPlot":
        model: Plot[Dataset, BoxPlotConfig] = Plot.initialize(
            plot_type=PlotType.BOX,
            pconfig=pconfig,
            anchor=anchor,
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
