import logging
from typing import Dict, List, Union

from plotly import graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType

logger = logging.getLogger(__name__)


# {'color': 'rgb(211,211,211,0.05)', 'name': 'background: EUR', 'x': -0.294, 'y': -1.527}
ElementT = Dict[str, Union[str, float, int]]


def plot(datasets: List[List[ElementT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    p = ScatterPlot(pconfig, len(datasets))

    return p.add_to_report(datasets, report)


class ScatterPlot(Plot):
    def __init__(self, pconfig: Dict, *args):
        super().__init__(PlotType.SCATTER, pconfig, *args)

        self.categories = pconfig.get("categories", [])
        self.marker_shape = {
            "size": 10,
            "line": {"width": 1},
            "opacity": 1,
        }

    def serialise(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().serialise()
        d["categories"] = self.categories
        d["marker_shape"] = self.marker_shape
        return d

    def populate_figure(
        self,
        fig: go.Figure,
        dataset: List[ElementT],
        is_log=False,
        is_pct=False,
        # is_flat=True,
    ) -> go.Figure:
        """
        Add traces to the figure
        """
        # Keeping track of element names, so we list them in legend only once
        met_names = set()

        for element in dataset:
            x = element["x"]

            if self.categories:
                if isinstance(x, int) and (0 <= x < len(self.categories)):
                    x = self.categories[x]
                else:
                    logger.error(
                        f"Scatter plot {self.id}: x={x} must be an index in the list of categories: {self.categories}"
                    )
                    continue

            marker = self.marker_shape.copy()
            if "marker_size" in element:
                marker["size"] = element["marker_size"]
            if "marker_line_width" in element:
                marker["line"]["width"] = element["marker_line_width"]
            if "color" in element:
                marker["color"] = element["color"]
            if "opacity" in element:
                marker["opacity"] = element["opacity"]

            show_in_legend = element["name"] not in met_names
            met_names.add(element["name"])

            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[element["y"]],
                    text=element["name"],
                    name=element["name"],
                    hoverinfo="text",
                    mode="markers",
                    marker=marker,
                    showlegend=show_in_legend,
                )
            )
        return fig

    def save_data_file(self, dataset: List, uid: str) -> None:
        pass
