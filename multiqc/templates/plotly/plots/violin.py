import dataclasses
import logging
from typing import Dict, List, Union, Optional

import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset

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
    data_by_sample: Dict[str, MetricsT]


class ViolinPlot(Plot):
    def __init__(
        self,
        pconfig: Dict,
        list_of_data_by_sample: List[Dict[str, MetricsT]],
        list_of_headers: List[Dict],
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_data_by_sample))

        # Extend each dataset object with a list of samples
        self.datasets: List[Dataset] = [
            Dataset(**d.__dict__, data_by_sample=self.prep_data(data_by_sample, headers))
            for d, data_by_sample, headers in zip(self.datasets, list_of_data_by_sample, list_of_headers)
        ]

        self.categories: List[str] = pconfig.get("categories", [])

        self.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            points="all",
            jitter=0,
        )

    @staticmethod
    def tt_label() -> Optional[str]:
        """Default tooltip label"""
        return "%{x}: %{y}"

    @staticmethod
    def prep_data(data_by_sample: Dict[str, MetricsT], headers: Dict) -> Dict[str, MetricsT]:
        # categories = []
        # for k, header in headers.items():
        #     categories.append(
        #         {
        #             "namespace": header["namespace"],
        #             "title": header["title"],
        #             "description": header["description"],
        #             "max": header["dmax"],
        #             "min": header["dmin"],
        #             "suffix": header.get("suffix", ""),
        #             "decimalPlaces": header.get("decimalPlaces", "2"),
        #         }
        #     )
        #     # Add the data
        #     thisdata = []
        #     these_snames = []
        #     for s_name, samp in data_by_sample.items():
        #         if k in samp:
        #             val = samp[k]
        #             dt.raw_vals[s_name][k] = val
        #
        #             if "modify" in header and callable(header["modify"]):
        #                 val = header["modify"](val)
        #
        #             thisdata.append(val)
        #             these_snames.append(s_name)

        # Reorder data to make it indexed by metric first, so each violin corresponds to
        # one metric, and contains a dot for each sample that have a value for this metric.
        # data_by_metric: Dict[str, MetricsT] = {}
        # for sample, metrics in data_by_sample.items():
        #     for metric, value in metrics.items():
        #         if metric not in data_by_metric:
        #             data_by_metric[metric] = {}
        #         data_by_metric[metric][sample] = value

        # replace metric key with header["title"]
        prep_metric_by_sample = dict()
        for sample, metrics_d in data_by_sample.items():
            prep_metric_by_sample[sample] = dict()
            for metric, value in metrics_d.items():
                if metric in headers:
                    title = headers[metric]["title"]
                    prep_metric_by_sample[sample][title] = value

        return prep_metric_by_sample

    def dump_for_javascript(self):
        """Serialise the data to pick up in plotly-js"""
        d = super().dump_for_javascript()
        d["categories"] = self.categories
        return d

    def create_figure(self, layout: go.Layout, dataset: Dataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)

        for name, metrics in dataset.data_by_sample.items():
            fig.add_trace(
                go.Violin(
                    name=name,
                    x=list(metrics.values()),
                    y=list(metrics.keys()),
                    **self.trace_params,
                )
            )
        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass
