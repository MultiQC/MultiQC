import logging
from typing import Dict, List, Union

from plotly import graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot

logger = logging.getLogger(__name__)


# {'color': 'rgb(211,211,211,0.05)', 'name': 'background: EUR', 'x': -0.294, 'y': -1.527}
ElementT = Dict[str, Union[str, float, int]]


def plot(
    datasets: List[List[ElementT]],
    pconfig: Dict,
) -> str:
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
        super().__init__("scatter", pconfig, *args)

    def populate_figure(self, fig: go.Figure, dataset: List[ElementT], is_log=False, is_pct=False) -> go.Figure:
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

            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[element["y"]],
                    text=element["name"],
                    hoverinfo="text",
                    mode="markers",
                    marker=dict(
                        size=element.get("marker_size", 20),
                        line_width=element.get("marker_line_width", 1),
                        color=element.get("color", "rgba(124, 181, 236, .5)"),
                        opacity=element.get("opacity", 1),
                    ),
                )
            )
        return fig

    def save_data_file(self, dataset: List, uid: str) -> None:
        pass

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()
        layout.update(
            {
                # "marker": {"color": "rgba(124, 181, 236, .5)"},
            }
        )
        return layout
