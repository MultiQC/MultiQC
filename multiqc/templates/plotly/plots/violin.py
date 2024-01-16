import dataclasses
import logging
from typing import Dict, List, Union, Any, Optional
import copy

import math
import numpy as np
import plotly.graph_objects as go

from multiqc.plots.table import make_table
from multiqc.plots.table_object import DataTable
from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset

logger = logging.getLogger(__name__)


def plot(dt: DataTable, show_table_by_default=False) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    p = ViolinPlot.from_dt(dt)

    from multiqc.utils import report

    plot_html = p.add_to_report(report)

    if show_table_by_default:
        table_html = make_table(dt, violin_switch=True)
        report.plot_data[p.id]["table_html"] = table_html
        report.plot_data[p.id]["plot_html"] = plot_html
        report.plot_data[p.id]["static"] = True
        return table_html
    else:
        return plot_html


THRESHOLD_BEFORE_OUTLIERS = 50


class ViolinPlot(Plot):
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
        ) -> "ViolinPlot.Dataset":
            outlier_indices_by_metric: Dict[str, List[int]] = dict()
            all_samples = set()
            show_only_outliers = False
            for i, metric in enumerate(headers_by_metric):
                header = headers_by_metric[metric]
                samples = samples_by_metric[metric]
                all_samples.update(set(samples))

                values = [v for v in values_by_metric[metric] if v is not None]
                if not values:
                    logger.warning(f"No non-None values for metric: {header['title']}")
                    values_by_metric[metric] = values
                    continue

                values_are_numerical = all(isinstance(v, (int, float)) for v in values)
                if values_are_numerical:
                    values = [v for v in values if not math.isnan(v)]
                    if not values:
                        logger.warning(f"All values are NaN for metric: {header['title']}")
                        values_by_metric[metric] = values
                        continue

                    xmin = header.get("min")
                    xmax = header.get("max")
                    if all(isinstance(v, (int, float)) for v in values):
                        if xmin is None or not isinstance(xmin, (int, float)):
                            xmin = min(values)
                        if xmax is None or not isinstance(xmin, (int, float)):
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

                        outlier_indices = find_outliers(
                            values,
                            minval=header.get("min"),
                            maxval=header.get("max"),
                        )

                        outlier_indices_by_metric[metric] = outlier_indices
                        logger.debug(f"Found {len(outlier_indices)} outliers for metric: {header['title']}")
                        show_only_outliers = True
                values_by_metric[metric] = values

            return ViolinPlot.Dataset(
                **dataset.__dict__,
                values_by_metric=values_by_metric,
                samples_by_metric=samples_by_metric,
                headers_by_metric=headers_by_metric,
                outlier_indices_by_metric=outlier_indices_by_metric,
                all_samples=list(all_samples),
                show_only_outliers=show_only_outliers,
            )

    def __init__(
        self,
        list_of_values_by_metric: List[Dict[str, List[Union[List[int], List[float], List[str]]]]],
        list_of_samples_by_metric: List[Dict[str, List[str]]],
        list_of_headers_by_metric: List[Dict[str, Dict]],
        pconfig: Dict,
        table_html: Optional[str] = None,
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_values_by_metric))

        self.table_html = table_html

        self.datasets: List[ViolinPlot.Dataset] = [
            ViolinPlot.Dataset.create(ds, values_by_metric, samples_by_metric, headers_by_metric)
            for ds, values_by_metric, samples_by_metric, headers_by_metric in zip(
                self.datasets,
                list_of_values_by_metric,
                list_of_samples_by_metric,
                list_of_headers_by_metric,
            )
        ]

        self.show_only_outliers = any(ds.show_only_outliers for ds in self.datasets)

        self.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            fillcolor="rgba(0,0,0,0.1)",
            points=False,  # Don't show points, we'll add them manually
            # The hover information is useful, but the formatting is ugly and not
            # configurable as far as I can see. Also, it's not possible to disable it,
            # so setting it to "points" as we don't show points, so it's effectively
            # disabling it.
            hoveron="points",
        )

        self.scatter_trace_params = {
            "mode": "markers",
            "marker": {"size": 4, "color": "rgba(0,0,0,1)"},
            "showlegend": False,
            "hovertemplate": self.trace_params["hovertemplate"],
            "hoverlabel": {"bgcolor": "white"},
        }

        num_rows = max(len(ds.values_by_metric) for ds in self.datasets)

        self.layout.update(
            height=70 * max(len(ds.values_by_metric) for ds in self.datasets),
            margin=dict(
                pad=0,
                # t=43,
                b=40,
            ),
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
    def from_dt(dt: DataTable, show_table_by_default: bool = False) -> "ViolinPlot":
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

        return ViolinPlot(
            [values_by_metric],
            [samples_by_metric],
            [headers_by_metric],
            pconfig=dt.pconfig,
        )

    @staticmethod
    def tt_label() -> str:
        return ": %{x}"

    # def interactive_plot(self, report) -> str:
    #     plot_html = super().interactive_plot(report)
    #     if self.table_html:
    #         self.plot_html = plot_html
    #         return self.table_html
    #     else:
    #         self.table_html = plot_html
    #         super().interactive_plot(report)

    def dump_for_javascript(self):
        """Serialise the data to pick up in plotly-js"""
        d = super().dump_for_javascript()
        d.update(
            {
                "scatter_trace_params": self.scatter_trace_params,
                "show_only_outliers": any(ds.show_only_outliers for ds in self.datasets),
                "initial_html": self.table_html,
            }
        )
        return d

    # def control_panel(self):
    #     """Add a control panel to the plot"""
    #     if not self.show_only_outliers:
    #         return ""
    #     return (
    #         f'<div class="beeswarm-hovertext" id="{self.id}-hoverinfo">' +
    #         '<em class="placeholder">Hover over a data point for more information</em>' +
    #         '</div>'
    #     )

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
                        xaxis=f"x{i + 1}",
                        yaxis=f"y{i + 1}",
                        showlegend=False,
                        **self.scatter_trace_params,
                    ),
                )
        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass


def find_outliers(
    values: Union[List[int], List[float]],
    n: Optional[int] = None,
    z_cutoff: float = 2.0,
    minval: Optional[Union[float, int]] = None,
    maxval: Optional[Union[float, int]] = None,
) -> List[int]:
    """
    If `n` is defined, find `n` most outlying points in a list.
    Otherwise, find outliers with a Z-score above `z_cutoff`.

    Return indices in the input list.

    `minval` and/or `maxval` can be added to include them into the outlier detection array.
    This is a trick to avoid the following problem: for example, percentage values are
    clustered around 0, but differ in e.g. 3+rd decimal place. In this case, the outliers
    will be nearly no different from the other values. If we add "100%" as an artificial
    additional value, none of those near-zero clustered values won't be called outliers.
    """
    if len(values) == 0 or (n is not None and n <= 0):
        return []

    if minval is not None:
        values = [minval] + values
    if maxval is not None:
        values = [maxval] + values
    values = np.array(values)

    # Calculate the mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    if std_dev == 0:
        logger.warning(f"All {len(values)} points have the same values")
        return []

    # Calculate Z-scores (measures of "outlyingness")
    z_scores = np.abs((values - mean) / std_dev)

    # Get indices of the top N outliers
    if n:
        indices = np.argsort(z_scores)[-n:].tolist()
    else:
        indices = []
        while z_cutoff > 1.0:
            indices = np.where(z_scores > z_cutoff)[0]
            if len(indices) > (int(minval is not None) + int(maxval is not None)):
                indices = indices.tolist()
                break
            logger.warning(f"No outliers found with Z-score cutoff {z_cutoff}, trying a lower threshold")
            z_cutoff -= 0.2
            indices = indices.tolist()

    if minval is not None and maxval is not None:
        if 0 in indices:
            indices.remove(0)
        if 1 in indices:
            indices.remove(1)
    elif minval is not None or maxval is not None:
        if 0 in indices:
            indices.remove(0)

    if len(indices) == 0:
        logger.warning(f"No outliers found with Z-score cutoff {z_cutoff}")
        indices = list(range(len(values)))

    return indices
