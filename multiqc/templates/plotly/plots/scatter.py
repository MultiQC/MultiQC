import copy
import dataclasses
import logging
from collections import defaultdict
from typing import Dict, List, Union, Optional

import numpy as np
from plotly import graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset

logger = logging.getLogger(__name__)


# {'color': 'rgb(211,211,211,0.05)', 'name': 'background: EUR', 'x': -0.294, 'y': -1.527}
PointT = Dict[str, Union[str, float, int]]


def plot(points_lists: List[List[PointT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param points_lists: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = ScatterPlot(pconfig, points_lists)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class Dataset(BaseDataset):
    points: List[PointT]


class ScatterPlot(Plot):
    def __init__(self, pconfig: Dict, points_lists: List[List[PointT]]):
        super().__init__(PlotType.SCATTER, pconfig, len(points_lists))

        # Extend each dataset object with a list of samples
        self.datasets: List[Dataset] = [
            Dataset(**d.__dict__, points=points) for d, points in zip(self.datasets, points_lists)
        ]

        self.categories = pconfig.get("categories", [])

        self.trace_params.update(
            textfont=dict(size=8),
            marker=dict(
                size=10,
                line=dict(width=1),
                opacity=1,
                color="rgba(124, 181, 236, .5)",
            ),
        )

        self.layout.height = self.layout.height or 600
        # Make a tooltip always show on hover over nearest point on plot
        self.layout.hoverdistance = -1

    @staticmethod
    def tt_label() -> Optional[str]:
        """Default tooltip label"""
        return "<br><b>X</b>: %{x}<br><b>Y</b>: %{y}>"

    def dump_for_javascript(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().dump_for_javascript()
        d["categories"] = self.categories
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
        fig = go.Figure(layout=layout)

        MAX_ANNOTATIONS = 10  # Maximum number of dots to be annotated directly on the plot
        n_annotated = len([el for el in dataset.points if "annotation" in el])
        if n_annotated < MAX_ANNOTATIONS:
            # Finding and marking outliers to only label them
            # 1. Calculate Z-scores
            points = [(i, x) for i, x in enumerate(dataset.points) if x.get("annotate", True) is not False]
            x_values = np.array([x["x"] for (i, x) in points])
            y_values = np.array([x["y"] for (i, x) in points])
            x_std = np.std(x_values)
            y_std = np.std(y_values)
            if x_std == 0 and y_std == 0:
                logger.warning(f"Scatter plot {self.id}: all {len(points)} points have the same coordinates")
                if len(points) == 1:  # Only single point - annotate it!
                    for i, point in points:
                        point["annotation"] = point["name"]
            else:
                x_z_scores = np.abs((x_values - np.mean(x_values)) / x_std) if x_std else np.zeros_like(x_values)
                y_z_scores = np.abs((y_values - np.mean(y_values)) / y_std) if y_std else np.zeros_like(y_values)
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
            n_annotated = len([point for point in dataset.points if "annotation" in point])

        # If there are few unique colors, we can additionally put a unique list into a legend
        # (even though some color might belong to many distinct names - we will just crop the list)
        names_by_legend_key = defaultdict(set)
        for el in dataset.points:
            legend_key = (el.get("color"), el.get("marker_size"), el.get("marker_line_width"), el.get("group"))
            names_by_legend_key[legend_key].add(el["name"])
        layout.showlegend = True

        in_legend = set()
        for element in dataset.points:
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
                    MAX_WIDTH = 60
                    if len(label) > MAX_WIDTH:  # Crop long labels
                        label = label[:MAX_WIDTH] + "..."
                    show_in_legend = True
                    name = label

            params = copy.deepcopy(self.trace_params)
            marker = params.pop("marker")
            if color:
                marker["color"] = color
            if "marker_line_width" in element:
                marker["line"]["width"] = element["marker_line_width"]
            if "marker_size" in element:
                marker["size"] = element["marker_size"]
            if "opacity" in element:
                marker["opacity"] = element["opacity"]

            if annotation:
                params["mode"] = "markers+text"
            if n_annotated > 0:  # Reduce opacity of the borders that clutter the annotations:
                marker["line"]["color"] = "rgba(0, 0, 0, .2)"

            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[element["y"]],
                    name=name,
                    text=[annotation or name],
                    showlegend=show_in_legend,
                    marker=marker,
                    **params,
                )
            )
        fig.layout.height += len(in_legend) * 5  # extra space for legend
        return fig

    def save_data_file(self, dataset: BaseDataset) -> None:
        pass