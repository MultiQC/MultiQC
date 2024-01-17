import dataclasses
import logging
from typing import Dict, List, Union, Any, Optional
import copy

import math
import numpy as np
import plotly.graph_objects as go

from multiqc.plots.table_object import DataTable
from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
from .table import make_table

logger = logging.getLogger(__name__)


def plot(dt: DataTable, show_table_by_default=False) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    p = ViolinPlot.from_dt(dt)

    from multiqc.utils import report

    violin_html = p.add_to_report(report)
    table_html, configuration_modal = make_table(dt, violin_switch=True)
    html = "".join(
        [
            f"<div id='mqc-violin-{p.id}' style='{'display: none;' if show_table_by_default else ''}'>{violin_html}</div>",
            f"<div id='mqc-table-{p.id}' style='{'display: none;' if not show_table_by_default else ''}'>{table_html}</div>",
            configuration_modal,
        ]
    )
    report.plot_data[p.id]["static"] = show_table_by_default

    return html


THRESHOLD_BEFORE_OUTLIERS = 50


class ViolinPlot(Plot):
    @dataclasses.dataclass
    class Dataset(BaseDataset):
        metrics: List[str]
        header_by_metric: Dict[str, Dict[str, Any]]
        values_by_sample_by_metric: Dict[str, Dict[str, Union[List[int], List[float], List[str]]]]
        outliers_by_metric: Dict[str, List[str]]
        all_samples: List[str]  # unique list of all samples in this dataset
        show_only_outliers: False

        @staticmethod
        def create(
            dataset: BaseDataset,
            values_by_sample_by_metric: Dict[str, Dict[str, Union[List[int], List[float], List[str]]]],
            headers_by_metric: Dict[str, Dict[str, Any]],
        ) -> "ViolinPlot.Dataset":
            ds = ViolinPlot.Dataset(
                **dataset.__dict__,
                metrics=list(headers_by_metric.keys()),
                header_by_metric=headers_by_metric,
                values_by_sample_by_metric=dict(),
                outliers_by_metric=dict(),
                all_samples=[],
                show_only_outliers=False,
            )

            all_samples = set()
            for metric in ds.metrics:
                header = ds.header_by_metric[metric]

                value_by_sample = values_by_sample_by_metric[metric]
                value_by_sample = {s: v for s, v in value_by_sample.items() if v is not None}
                if not value_by_sample:
                    logger.warning(f"No non-None values for metric: {header['title']}")
                    ds.values_by_sample_by_metric[metric] = value_by_sample
                    continue

                values_are_numerical = all(isinstance(v, (int, float)) for v in value_by_sample.values())
                if values_are_numerical:
                    value_by_sample = {s: v for s, v in value_by_sample.items() if not math.isnan(v)}
                    if not value_by_sample:
                        logger.warning(f"All values are NaN for metric: {header['title']}")
                        ds.values_by_sample_by_metric[metric] = value_by_sample
                        continue

                    xmin = header.get("min")
                    xmax = header.get("max")
                    if all(isinstance(v, (int, float)) for v in value_by_sample.values()):
                        if xmin is None or not isinstance(xmin, (int, float)):
                            xmin = min(value_by_sample.values())
                        if xmax is None or not isinstance(xmin, (int, float)):
                            xmax = max(value_by_sample.values())
                        xmin -= (xmax - xmin) * 0.005
                        xmax += (xmax - xmin) * 0.005
                    header["xaxis"] = {"range": [xmin, xmax]}

                    if len(value_by_sample) > THRESHOLD_BEFORE_OUTLIERS:
                        logger.warning(
                            f"Violin plot with {len(value_by_sample)} > {THRESHOLD_BEFORE_OUTLIERS} samples. "
                            f"This may be too many to display clearly, so showing "
                            f"only outliers in each violin."
                        )

                        samples = list(value_by_sample.values())
                        values = list(value_by_sample.values())
                        outlier_indices = find_outliers(
                            values,
                            minval=header.get("min"),
                            maxval=header.get("max"),
                        )
                        logger.debug(f"Found {len(outlier_indices)} outliers for metric: {header['title']}")
                        ds.outliers_by_metric[metric] = [samples[idx] for idx in outlier_indices]
                        ds.show_only_outliers = True

                ds.values_by_sample_by_metric[metric] = value_by_sample
                all_samples.update(set(list(value_by_sample.keys())))

            ds.all_samples = sorted(all_samples)
            return ds

    def __init__(
        self,
        list_of_values_by_sample_by_metric: List[Dict[str, Dict[str, Union[List[int], List[float], List[str]]]]],
        list_of_header_by_metric: List[Dict[str, Dict]],
        pconfig: Dict,
        table_html: Optional[str] = None,
    ):
        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_values_by_sample_by_metric))

        self.table_html = table_html

        self.datasets: List[ViolinPlot.Dataset] = [
            ViolinPlot.Dataset.create(ds, values_by_sample_by_metric, headers_by_metric)
            for ds, values_by_sample_by_metric, headers_by_metric in zip(
                self.datasets,
                list_of_values_by_sample_by_metric,
                list_of_header_by_metric,
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

        total_rows = 0
        unhidden_rows = 0
        for ds in self.datasets:
            total_rows = max(total_rows, len(ds.metrics))
            unhidden_metrics = [m for m in ds.metrics if not ds.header_by_metric[m].get("hidden", False)]
            unhidden_rows = max(unhidden_rows, len(unhidden_metrics))

        self.layout.update(
            height=70 * unhidden_rows,
            margin=dict(
                pad=0,
                # t=43,
                b=40,
            ),
            violingap=0,
            grid=dict(
                rows=unhidden_rows,
                columns=1,
                roworder="top to bottom",
                ygap=0.4,
                subplots=[[(f"x{i + 1}y{i + 1}" if i > 0 else "xy")] for i in range(total_rows)],
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
    def from_dt(dt: DataTable) -> "ViolinPlot":
        values_by_sample_by_metric = dict()
        header_by_metric = dict()

        for idx, k, header in dt.get_headers_in_order():
            rid = header["rid"]
            header_by_metric[rid] = {
                "rid": header["rid"],
                "namespace": header["namespace"],
                "title": header["title"],
                "description": header["description"],
                "max": header.get("max"),
                "min": header.get("min"),
                "suffix": header.get("suffix", ""),
                "color": header.get("color"),
                "hidden": header.get("hidden"),
            }
            values_by_sample_by_metric[rid] = dict()
            for s_name, val_by_metric in dt.data[idx].items():
                if k in val_by_metric:
                    values_by_sample_by_metric[rid][s_name] = val_by_metric[k]

        # for idx, (_data_by_sample, _headers_by_metric) in enumerate(zip(dt.data, dt.headers)):
        #     for sample, data in _data_by_sample.items():
        #         for metric, value in data.items():
        #             if metric not in _headers_by_metric:
        #                 continue
        #             if metric not in list_of_values_by_sample_by_metric:
        #                 list_of_values_by_sample_by_metric[metric] = dict()
        #             list_of_values_by_sample_by_metric[metric][sample] = value
        #
        #     for metric, header in _headers_by_metric.items():
        #         color = header.get("colour")
        #         headers_by_metric[metric] = {
        #             "namespace": header["namespace"],
        #             "title": header["title"],
        #             "description": header["description"],
        #             "max": header.get("max"),
        #             "min": header.get("min"),
        #             "suffix": header.get("suffix", ""),
        #             "color": color,
        #         }

        return ViolinPlot(
            [values_by_sample_by_metric],
            [header_by_metric],
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

    def buttons(self) -> []:
        """Add a control panel to the plot"""
        buttons = [
            self._btn(
                cls="mqc-violin-to-table",
                label="<span class='glyphicon glyphicon-th-list'></span> Switch to table",
            )
        ]
        if any(len(ds.metrics) > 1 for ds in self.datasets):
            buttons.append(
                self._btn(
                    cls="mqc_table_configModal_btn",
                    label="<span class='glyphicon glyphicon-th'></span> Configure columns",
                    data_attrs={"toggle": "modal", "target": f"#{self.id}_configModal"},
                )
            )

        return buttons + super().buttons()

    def create_figure(self, layout: go.Layout, dataset: Dataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """
        for i, header in enumerate(dataset.header_by_metric):
            layout[f"xaxis{i + 1}"] = copy.deepcopy(layout["xaxis"])
            layout[f"xaxis{i + 1}"].update(header["xaxis"])
            layout[f"yaxis{i + 1}"] = copy.deepcopy(layout["yaxis"])

        layout.showlegend = False
        fig = go.Figure(layout=layout)

        for i, metric in enumerate(dataset.header_by_metric):
            header = dataset.header_by_metric[metric]
            values_by_sample = dataset.values_by_sample_by_metric[metric]
            params = copy.deepcopy(self.trace_params)
            if header.get("color"):
                params["fillcolor"] = f"rgba({header['color']},0.5)"
            fig.add_trace(
                go.Violin(
                    x=list(values_by_sample.values()),
                    name=header["title"] + "  ",
                    text=list(values_by_sample.keys()),
                    xaxis=f"x{i + 1}",
                    yaxis=f"y{i + 1}",
                    **params,
                ),
            )
            for j, (sample, value) in enumerate(values_by_sample.items()):
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
