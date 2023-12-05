import logging
from collections import defaultdict
from typing import Dict, List, Union

import numpy as np
from plotly import graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, Dataset

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
    p = ScatterPlot(pconfig, datasets)

    from multiqc.utils import report

    return p.add_to_report(report)


class ScatterPlot(Plot):
    def __init__(self, pconfig: Dict, datasets: List):
        super().__init__(PlotType.SCATTER, pconfig, datasets)

        self.height = pconfig.get("height", 600)
        self.categories = pconfig.get("categories", [])
        self.default_marker = {
            "size": 10,
            "line": {"width": 1},
            "opacity": 1,
            "color": "rgba(124, 181, 236, .5)",
        }

    def serialise(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().serialise()
        d["categories"] = self.categories
        d["default_marker"] = self.default_marker
        return d

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

        import json

        with open(f"/Users/vlad/git/playground/{self.id}-layout.json", "w") as f:
            f.write(json.dumps(layout.to_plotly_json()))
        with open(f"/Users/vlad/git/playground/{self.id}-data.json", "w") as f:
            f.write(json.dumps(dataset.data))

        data: List[ElementT] = dataset.data

        fig = go.Figure(layout=layout)

        # On flat plots, we try to annotate some points directly on the plot: outliers of if there are only few
        show_annotations = False

        names_by_color = defaultdict(set)
        for el in data:
            names_by_color[el.get("color")].add(el["name"])

        if len(data) <= 10:
            # Only few data points: we can annotate all elements directly on the plot.
            layout.showlegend = False  # No need to add a legend
            show_annotations = True
            for element in data:
                element["annotation"] = element["name"]
        else:
            # Too many data points, finding and marking outliers to only label them
            show_annotations = True

            # Finding outliers: calculate Z-scores
            x_values = np.array([point["x"] for point in data])
            y_values = np.array([point["y"] for point in data])
            x_z_scores = np.abs((x_values - np.mean(x_values)) / np.std(x_values))
            y_z_scores = np.abs((y_values - np.mean(y_values)) / np.std(y_values))

            # Somewhat arbitrary threshold for outliers
            OUTLIER_THRESHOLD = 2

            for i, element in enumerate(data):
                # Check if point is an outlier or if total points are less than 10
                if x_z_scores[i] > OUTLIER_THRESHOLD or y_z_scores[i] > OUTLIER_THRESHOLD:
                    element["annotation"] = element["name"]

            if len(names_by_color) <= 10:
                # There are a lot of dots, but few unique colors. We can additionally put
                # a unique list into a legend (even though some color might belong to many
                # distinct names - we will just crop the list)
                layout.showlegend = True

        in_legend = set()
        for element in data:
            x = element["x"]
            if self.categories:
                if isinstance(x, float):
                    x = int(round(x))
                if isinstance(x, int) and (0 <= x < len(self.categories)):
                    x = self.categories[x]
                else:
                    logger.error(
                        f"Scatter plot {self.id}: x={x} must be an index in the list of categories: {self.categories}"
                    )
                    continue

            name = element["name"]
            annotation = element.get("annotation")
            color = element.get("color")

            show_in_legend = False
            if layout.showlegend and not element.get("hide_in_legend"):
                if color not in in_legend:
                    in_legend.add(color)
                    show_in_legend = True
                    name = ", ".join(sorted(names_by_color[color]))
                    if len(name) > 50:  # Crop long names
                        name = name[:50] + "..."

            marker = self.default_marker.copy()
            # if not show_annotations:
            #     # Removing default color to make sure Plotly assigns distinct colors to distinct names
            #     del marker["color"]
            if color:
                marker["color"] = color
            if show_annotations:
                marker["line"]["width"] = 0  # Remove the borders that clutter the annotations
            elif "marker_line_width" in element:
                marker["line"]["width"] = element["marker_line_width"]
            if "marker_size" in element:
                marker["size"] = element["marker_size"]
            if "opacity" in element:
                marker["opacity"] = element["opacity"]

            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[element["y"]],
                    name=name,
                    text=annotation,
                    textfont=dict(size=8),
                    mode="markers+text" if annotation else "markers",
                    marker=marker,
                    showlegend=show_in_legend,
                )
            )
        return fig

    def save_data_file(self, dataset: Dataset) -> None:
        pass
