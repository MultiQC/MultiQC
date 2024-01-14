import dataclasses
import logging
from typing import Dict, List, Union, Any, Optional
import copy

import numpy as np
import plotly.graph_objects as go

from multiqc.plots.table_object import DataTable
from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset

logger = logging.getLogger(__name__)


def plot(dt: DataTable) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    values_by_metric = dict()
    samples_by_metric = dict()
    headers_by_metric = dict()

    for idx, (_data_by_sample, _headers_by_metric) in enumerate(zip(dt.data, dt.headers)):
        for sample, data in _data_by_sample.items():
            for metric, value in data.items():
                if metric not in _headers_by_metric:
                    continue
                if metric not in values_by_metric:
                    values_by_metric[metric] = []
                if metric not in samples_by_metric:
                    samples_by_metric[metric] = []
                samples_by_metric[metric].append(sample)
                values_by_metric[metric].append(value)

        for metric, header in _headers_by_metric.items():
            color = header.get("colour")
            headers_by_metric[metric] = {
                "namespace": header["namespace"],
                "title": header["title"],
                "description": header["description"],
                "max": header.get("max"),
                "min": header.get("min"),
                "suffix": header.get("suffix", ""),
                "color": color,
            }

    p = ViolinPlot(
        [values_by_metric],
        [samples_by_metric],
        [headers_by_metric],
        dt.pconfig,
    )

    from multiqc.utils import report

    return p.add_to_report(report)


THRESHOLD_BEFORE_OUTLIERS = 50
# NUMBER_OF_OUTLIERS = 20


@dataclasses.dataclass
class Dataset(BaseDataset):
    values_by_metric: Dict[str, Union[List[int], List[float], List[str]]]
    samples_by_metric: Dict[str, List[str]]
    headers_by_metric: Dict[str, Dict[str, Any]]
    outlier_indices_by_metric: Dict[str, List[int]]
    all_samples: List[str]  # list of all samples in this dataset
    show_only_outliers: False

    @staticmethod
    def create(
        dataset: BaseDataset,
        values_by_metric: Dict[str, Union[List[int], List[float], List[str]]],
        samples_by_metric: Dict[str, List[str]],
        headers_by_metric: Dict[str, Dict[str, Any]],
    ) -> "Dataset":
        outlier_indices_by_metric: Dict[str, List[int]] = dict()
        all_samples = set()
        show_only_outliers = False
        for i, metric in enumerate(headers_by_metric):
            header = headers_by_metric[metric]
            values = values_by_metric[metric]
            samples = samples_by_metric[metric]

            all_samples.update(set(samples))

            xmin = header.get("min")
            if xmin is None:
                xmin = min(values)
            xmax = header.get("max")
            if xmax is None:
                xmax = max(values)
            xmin -= (xmax - xmin) * 0.005
            xmax += (xmax - xmin) * 0.005
            header["xaxis"] = {"range": [xmin, xmax]}

            if len(values) > THRESHOLD_BEFORE_OUTLIERS:
                logger.warning(
                    f"Violin plot with {len(values)} > {THRESHOLD_BEFORE_OUTLIERS} samples. "
                    f"This may be too many to display clearly, so showing "
                    f"only outliers in each violin."
                )
                outlier_indices = find_outliers(values)
                outlier_indices_by_metric[metric] = outlier_indices
                logger.debug(f"Found {len(outlier_indices)} outliers for metric: {header['title']}")
                show_only_outliers = True

        return Dataset(
            **dataset.__dict__,
            values_by_metric=values_by_metric,
            samples_by_metric=samples_by_metric,
            headers_by_metric=headers_by_metric,
            outlier_indices_by_metric=outlier_indices_by_metric,
            all_samples=list(all_samples),
            show_only_outliers=show_only_outliers,
        )


class ViolinPlot(Plot):
    def __init__(
        self,
        list_of_values_by_metric: List[Dict[str, List[Union[List[int], List[float], List[str]]]]],
        list_of_samples_by_metric: List[Dict[str, List[str]]],
        list_of_headers_by_metric: List[Dict[str, Dict]],
        pconfig: Dict,
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_values_by_metric))

        self.datasets: List[Dataset] = [
            Dataset.create(ds, values_by_metric, samples_by_metric, headers_by_metric)
            for ds, values_by_metric, samples_by_metric, headers_by_metric in zip(
                self.datasets,
                list_of_values_by_metric,
                list_of_samples_by_metric,
                list_of_headers_by_metric,
            )
        ]

        self.categories: List[str] = pconfig.get("categories", [])

        self.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            fillcolor="rgba(0,0,0,0.1)",
            line={"width": 0},  # Disable the border (still have the fill color)
            points=False,  # Don't show points, we'll add them manually
            # The hover information is useful, but the formatting is ugly and not
            # configurable as far as I can see. Also, it's not possible to disable it,
            # so setting it to "points" as we don't show points, so it's effectively
            # disabling it.
            hoveron="points",
        )

        num_rows = max(len(ds.values_by_metric) for ds in self.datasets)

        self.layout.update(
            height=70 * max(len(ds.values_by_metric) for ds in self.datasets),
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

    def add_to_report(self, report) -> str:
        show_only_outliers = any(ds.show_only_outliers for ds in self.datasets)
        warning = ""
        if show_only_outliers:
            warning = (
                f'<p class="text-muted">The number of points is above {THRESHOLD_BEFORE_OUTLIERS}, '
                f"so for efficiency, separate points are shown only for outliers.</p>"
            )
        return warning + super().add_to_report(report)

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
        for i, header in enumerate(dataset.headers_by_metric):
            layout[f"xaxis{i + 1}"] = copy.deepcopy(layout["xaxis"])
            layout[f"xaxis{i + 1}"].update(header["xaxis"])
            layout[f"yaxis{i + 1}"] = copy.deepcopy(layout["yaxis"])

        layout.showlegend = False
        fig = go.Figure(layout=layout)

        for i, metric in enumerate(dataset.headers_by_metric):
            header = dataset.headers_by_metric[metric]
            values = dataset.values_by_metric[metric]
            samples = dataset.samples_by_metric[metric]
            params = copy.deepcopy(self.trace_params)
            if header.get("color"):
                params["fillcolor"] = f"rgba({header['color']},0.5)"
            fig.add_trace(
                go.Violin(
                    x=values,
                    name=header["title"] + "  ",
                    text=samples,
                    xaxis=f"x{i + 1}",
                    yaxis=f"y{i + 1}",
                    **params,
                ),
            )
            for j, (sample, value) in enumerate(zip(samples, values)):
                fig.add_trace(
                    go.Scatter(
                        x=[value],
                        y=[header["title"] + "  "],
                        text=[sample],
                        mode="markers",
                        xaxis=f"x{i + 1}",
                        yaxis=f"y{i + 1}",
                        showlegend=False,
                        hovertemplate=self.trace_params["hovertemplate"],
                    ),
                )
        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass


def find_outliers(values: Union[List[int], List[float]], n: Optional[int] = None, z_cutoff: float = 3.0) -> List[int]:
    """
    If `n` is defined, find `n` most outlying points in a list.
    Otherwise, find outliers with a Z-score above `z_cutoff`.

    Return indices in the input list.
    """
    if len(values) == 0 or (n is not None and n <= 0):
        return []

    values = np.array(values)

    # Calculate the mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    if std_dev == 0:
        logger.warning(f"All {len(values)} points have the same values, just returning the first point")
        return [0]

    # Calculate Z-scores (measures of "outlyingness")
    z_scores = np.abs((values - mean) / std_dev)

    # Get indices of the top N outliers
    if n:
        outliers_indices = np.argsort(z_scores)[-n:]
        return outliers_indices.tolist()
    else:
        return np.where(z_scores > z_cutoff)[0].tolist()
