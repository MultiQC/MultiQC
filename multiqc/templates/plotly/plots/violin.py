import dataclasses
import logging
from typing import Dict, List, Union
import copy

import plotly.graph_objects as go

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
    samples: List[str]  # list of all samples in this dataset
    sample_colors: Dict[str, str]  # a color matching each sample


class ViolinPlot(Plot):
    def __init__(
        self,
        pconfig: Dict,
        list_of_data_by_sample: List[Dict[str, MetricsT]],
        list_of_header_by_metric: List[Dict],
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_data_by_sample))

        c_scale = mqc_colour.mqc_colour_scale("plot_defaults")

        # Extend each dataset object with a list of samples
        datasets: List[Dataset] = []
        for ds, data_by_sample, header_by_metric in zip(
            self.datasets, list_of_data_by_sample, list_of_header_by_metric
        ):
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

            for i, (metric, header) in enumerate(header_by_metric.items()):
                data = data_by_metric[metric]
                xmin = header.get("min")
                if xmin is None:
                    xmin = min(data.values())
                    xmin -= xmin * 0.05
                xmax = header_by_metric[metric].get("max")
                if xmax is None:
                    xmax = max(data.values())
                    xmax += xmax * 0.05
                header_by_metric[metric]["xaxis"] = {
                    "rangemode": "tozero" if xmin == 0 else "normal",
                    "range": [xmin, xmax],
                }
                # header_by_metric[metric]["color"] = "rgba(0,0,0,0.5)"

            all_samples = list(data_by_sample.keys())
            sample_colors = {sn: c_scale.get_colour(i, lighten=1) for i, sn in enumerate(all_samples)}
            datasets.append(
                Dataset(
                    *ds.__dict__,
                    data_by_metric=data_by_metric,
                    header_by_metric=header_by_metric,
                    samples=all_samples,
                    sample_colors=sample_colors,
                )
            )

        self.datasets = datasets

        self.categories: List[str] = pconfig.get("categories", [])

        self.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            # jitter=0.5,
            # pointpos=0,
            # points="all",
            fillcolor="rgba(0,0,0,0.1)",
            line={"width": 0},
            marker={"color": "rgb(55,126,184)", "size": 4},
            # The hover information is useful, but the formatting is ugly and not
            # configurable as far as I can see. Also, it's not possible to disable it,
            # so setting it to "points" as we don't show points, so it's effectively
            # disabling it.
            hoveron="points",
        )

        num_rows = max(len(ds.header_by_metric) for ds in self.datasets)

        self.layout.update(
            height=70 * max(len(ds.header_by_metric) for ds in self.datasets),
            margin=dict(pad=0, t=10, b=30),
            violingap=0,
            grid=dict(
                rows=num_rows,
                columns=1,
                roworder="top to bottom",
                ygap=0.4,
                subplots=[[(f"x{i + 1}y{i + 1}" if i > 0 else "xy")] for i in range(num_rows)],
            ),
            xaxis=dict(
                tickfont=dict(size=9, color="rgba(0,0,0,0.5)"),
            ),
            yaxis=dict(
                automargin=True,
            ),
            # hoverlabel=dict(
            #     bgcolor="rgb(141,203,255)",
            # ),
        )

    @staticmethod
    def tt_label() -> str:
        return ": %{x}"

    def dump_for_javascript(self):
        """Serialise the data to pick up in plotly-js"""
        d = super().dump_for_javascript()
        d["categories"] = self.categories
        return d

    def create_figure(self, layout: go.Layout, dataset: Dataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """
        for i, header in enumerate(dataset.header_by_metric.values()):
            layout[f"xaxis{i + 1}"] = copy.deepcopy(layout["xaxis"])
            layout[f"xaxis{i + 1}"].update(header["xaxis"])
            layout[f"yaxis{i + 1}"] = copy.deepcopy(layout["yaxis"])

        layout.showlegend = False
        fig = go.Figure(layout=layout)

        for i, (metric, header) in enumerate(dataset.header_by_metric.items()):
            data = dataset.data_by_metric[metric]
            fig.add_trace(
                go.Violin(
                    x=list(data.values()),
                    name=header.get("title", metric) + "  ",
                    text=list(data.keys()),
                    xaxis=f"x{i + 1}",
                    yaxis=f"y{i + 1}",
                    **self.trace_params,
                ),
            )
            for j, (sample, value) in enumerate(data.items()):
                fig.add_trace(
                    go.Scatter(
                        x=[value],
                        y=[header.get("title", metric) + "  "],
                        text=[sample],
                        mode="markers",
                        marker=dict(
                            color=dataset.sample_colors[sample],
                        ),
                        xaxis=f"x{i + 1}",
                        yaxis=f"y{i + 1}",
                        showlegend=False,
                        hovertemplate=self.trace_params["hovertemplate"],
                    ),
                )
        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass
