import logging
from typing import Dict, List, Union
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, Dataset

logger = logging.getLogger(__name__)


ElemT = Union[str, float, int]


def plot(
    dataset: List[List[ElemT]],
    xcats: List[str],
    ycats: List[str],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param dataset: One dataset. A dataset is a list of rows of values
    :param xcats: Labels for X axis
    :param ycats: Labels for Y axis
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = HeatmapPlot(pconfig, [dataset], xcats, ycats)

    from multiqc.utils import report

    return p.add_to_report(report)


class HeatmapPlot(Plot):
    def __init__(self, pconfig: Dict, datasets: List, xcats: List[str], ycats: List[str]):
        super().__init__(PlotType.HEATMAP, pconfig, datasets)
        self.xcats = xcats
        self.ycats = ycats

        if not self.height:
            self.height = len(self.ycats) * 17
            self.height = max(500, self.height)
        if not self.width:
            self.width = self.height

    def serialise(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().serialise()
        d["xcats"] = self.xcats
        d["ycats"] = self.ycats
        return d

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()
        layout.update(
            {
                "xaxis_tickangle": 45,
                "xaxis_nticks": len(self.xcats),
                "yaxis_nticks": len(self.ycats),
            }
        )
        return layout

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
        data: List[List[ElemT]] = dataset.data
        fig = go.Figure(
            data=go.Heatmap(
                z=data,
                x=self.xcats,
                y=self.ycats,
            ),
            layout=layout,
        )
        return fig

    def save_data_file(self, dataset: Dataset) -> None:
        pass
