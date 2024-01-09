import dataclasses
import logging
from typing import Dict, List, Union, Optional

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
from multiqc.utils import mqc_colour

logger = logging.getLogger(__name__)


# {"metric1": 0.1, "metric2": "v2"...}
MetricsT = Dict[str, Union[float, int, str]]


def plot(data: List[Dict[str, MetricsT]], headers: List[Dict], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    p = ViolinPlot(pconfig, data, headers)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class Dataset(BaseDataset):
    data_by_metric: Dict[str, MetricsT]
    header_by_metric: Dict[str, Dict[str, Union[str, int, float]]]


class ViolinPlot(Plot):
    def __init__(
        self,
        pconfig: Dict,
        list_of_data_by_sample: List[Dict[str, MetricsT]],
        list_of_header_by_metric: List[Dict],
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_data_by_sample))

        # Extend each dataset object with a list of samples
        for idx, (data_by_sample, header_by_metric) in enumerate(zip(list_of_data_by_sample, list_of_header_by_metric)):
            header_by_metric = {
                m: {k: v for k, v in header.items() if isinstance(v, (str, int, float))}
                for m, header in header_by_metric.items()
            }

            data_by_metric = {}
            for sample, metrics in data_by_sample.items():
                for metric, value in list(metrics.items()):
                    if metric not in header_by_metric:
                        continue
                    if metric not in data_by_metric:
                        data_by_metric[metric] = {}
                    data_by_metric[metric][sample] = value

            for metric, data in data_by_metric.items():
                xmin = header_by_metric[metric].get("min")
                header_by_metric[metric]["tozero"] = xmin == 0
                if xmin is None:
                    xmin = min(data.values())
                    xmin -= xmin * 0.05
                header_by_metric[metric]["min"] = xmin
                xmax = header_by_metric[metric].get("max")
                if xmax is None:
                    xmax = max(data.values())
                    xmax += xmax * 0.05
                header_by_metric[metric]["max"] = xmax

            self.datasets[idx] = Dataset(
                *self.datasets[idx].__dict__, data_by_metric=data_by_metric, header_by_metric=header_by_metric
            )

        self.categories: List[str] = pconfig.get("categories", [])

        self.trace_params.update(
            orientation="h",
            box={"visible": False},
            meanline={"visible": True},
            jitter=0.5,
            points="all",
            pointpos=0,
        )

        self.layout.height = 70 * max(len(ds.header_by_metric) for ds in self.datasets)
        self.layout.margin = dict(pad=0)
        self.layout.violingap = 0

    def dump_for_javascript(self):
        """Serialise the data to pick up in plotly-js"""
        d = super().dump_for_javascript()
        d["categories"] = self.categories
        return d

    @staticmethod
    def tt_label() -> Optional[str]:
        """Default tooltip label"""
        return "%{x}: %{y}"

    def create_figure(self, layout: go.Layout, dataset: Dataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """
        c_scale = mqc_colour.mqc_colour_scale("plot_defaults")
        fig = go.Figure(layout=layout)
        fig = make_subplots(
            len(dataset.data_by_metric),
            1,
            figure=fig,
            vertical_spacing=0.5 / len(dataset.data_by_metric),
        )
        for idx, (metric, data) in enumerate(dataset.data_by_metric.items()):
            fig.add_trace(
                go.Violin(
                    x=list(data.values()),
                    name=dataset.header_by_metric[metric].get("title", metric),
                    text=list(data.keys()),
                    marker={"color": c_scale.get_colour(idx, lighten=1)},
                    **self.trace_params,
                ),
                row=idx + 1,
                col=1,
            )
            fig.update_xaxes(
                range=[dataset.header_by_metric[metric]["min"], dataset.header_by_metric[metric]["max"]],
                rangemode="tozero" if dataset.header_by_metric[metric]["tozero"] else "normal",
                row=idx + 1,
                col=1,
            )

        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass
