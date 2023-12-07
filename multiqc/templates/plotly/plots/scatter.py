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

        self.layout.height = self.layout.height or 600
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

        # import json
        #
        # with open(f"/Users/vlad/git/playground/{self.id}-layout.json", "w") as f:
        #     f.write(json.dumps(layout.to_plotly_json()))
        # with open(f"/Users/vlad/git/playground/{self.id}-data.json", "w") as f:
        #     f.write(json.dumps(dataset.data))

        data: List[ElementT] = dataset.data

        fig = go.Figure(layout=layout)

        MAX_ANNOTATIONS = 10  # Maximum number of dots to be annotated directly on the plot
        n_annotated = len([el for el in data if "annotation" in el])
        if n_annotated < MAX_ANNOTATIONS:
            # Finding and marking outliers to only label them
            # 1. Calculate Z-scores
            points = [(i, x) for i, x in enumerate(data) if x.get("annotate", True) is not False]
            x_values = np.array([x["x"] for (i, x) in points])
            y_values = np.array([x["y"] for (i, x) in points])
            x_z_scores = np.abs((x_values - np.mean(x_values)) / np.std(x_values))
            y_z_scores = np.abs((y_values - np.mean(y_values)) / np.std(y_values))

            # 2. Find a reasonable threshold so there are not too many outliers
            threshold = 1.0
            while threshold <= 6.0:
                n_outliers = np.count_nonzero((x_z_scores > threshold) | (y_z_scores > threshold))
                logger.debug(f"Scatter plot outlier threshold: {threshold}, outliers: {n_outliers}")
                if n_annotated + n_outliers <= MAX_ANNOTATIONS:
                    break
                # If there are too many outliers, we increase the threshold until we have less than 10
                threshold += 0.2

            # 3. Annotate outliers that pass the threshold
            for (i, point), x_z_score, y_z_score in zip(points, x_z_scores, y_z_scores):
                # Check if point is an outlier or if total points are less than 10
                if x_z_score > threshold or y_z_score > threshold:
                    point["annotation"] = point["name"]
            n_annotated = len([point for point in data if "annotation" in point])

        # If there are few unique colors, we can additionally put a unique list into a legend
        # (even though some color might belong to many distinct names - we will just crop the list)
        names_by_legend_key = defaultdict(set)
        for el in data:
            legend_key = (el.get("color"), el.get("marker_size"), el.get("marker_line_width"), el.get("group"))
            names_by_legend_key[legend_key].add(el["name"])
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
            group = element.get("group")
            color = element.get("color")
            annotation = element.get("annotation")

            show_in_legend = False
            if layout.showlegend and not element.get("hide_in_legend"):
                key = (color, element.get("marker_size"), element.get("marker_line_width"), group)
                if key not in in_legend:
                    in_legend.add(key)
                    names = sorted(names_by_legend_key.get(key))
                    label = ", ".join(names)
                    if group:
                        label = f"{group}: {label}"
                    if len(label) > 70:  # Crop long labels
                        label = label[:70] + "..."
                    show_in_legend = True
                    name = label

            marker = self.default_marker.copy()
            if color:
                marker["color"] = color
            if n_annotated > 0:
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
