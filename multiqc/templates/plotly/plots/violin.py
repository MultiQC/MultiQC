import dataclasses
import logging
from typing import Dict, List, Union, Any, Optional
import copy

import math
import numpy as np
import plotly.graph_objects as go

from multiqc.plots.table_object import DataTable
from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
from multiqc.utils import config
from .table import make_table

logger = logging.getLogger(__name__)


def plot(dt: DataTable, show_table_by_default=False) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    p = ViolinPlot.from_dt(dt, show_table_by_default)

    from multiqc.utils import report

    return p.add_to_report(report)


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
            header_by_metric: Dict[str, Dict[str, Any]],
        ) -> "ViolinPlot.Dataset":
            ds = ViolinPlot.Dataset(
                **dataset.__dict__,
                metrics=list(header_by_metric.keys()),
                header_by_metric=header_by_metric,
                values_by_sample_by_metric=dict(),
                outliers_by_metric=dict(),
                all_samples=[],
                show_only_outliers=False,
            )

            all_samples = set()
            for metric in ds.metrics:
                header = ds.header_by_metric[metric]

                header["xaxis"] = {"ticksuffix": header.get("suffix")}

                value_by_sample = values_by_sample_by_metric[metric]
                value_by_sample = {s: v for s, v in value_by_sample.items() if v is not None}
                for s, v in value_by_sample.items():
                    if isinstance(v, str):
                        try:
                            value_by_sample[s] = int(v)
                        except ValueError:
                            try:
                                value_by_sample[s] = float(v)
                            except ValueError:
                                pass
                    if "modify" in header and callable(header["modify"]):
                        try:
                            value_by_sample[s] = header["modify"](v)
                        except Exception:
                            pass

                if not value_by_sample:
                    logger.warning(f"No non-None values for metric: {header['title']}")
                    ds.values_by_sample_by_metric[metric] = value_by_sample
                    continue

                values_are_numeric = all(isinstance(v, (int, float)) for v in value_by_sample.values())
                values_are_integer = all(isinstance(v, int) for v in value_by_sample.values())

                if values_are_numeric:
                    value_by_sample = {s: v for s, v in value_by_sample.items() if not math.isnan(v)}
                    if not value_by_sample:
                        logger.warning(f"All values are NaN for metric: {header['title']}")
                        ds.values_by_sample_by_metric[metric] = value_by_sample
                        continue

                    xmin = header.get("dmin")
                    xmax = header.get("dmax")
                    tickvals = None
                    if xmin == xmax == 0:  # Plotly will modify the 0:0 range to -1:1, and we want to keep it 0-centered
                        xmax = 1
                        tickvals = [0]
                    header["xaxis"]["tickvals"] = tickvals
                    if xmin is not None and xmax is not None:
                        xmin -= (xmax - xmin) * 0.005
                        xmax += (xmax - xmin) * 0.005
                        header["xaxis"]["range"] = [xmin, xmax]

                    if len(value_by_sample) > THRESHOLD_BEFORE_OUTLIERS:
                        logger.debug(
                            f"Violin for '{header['title']}': sample number is {len(value_by_sample)} > {THRESHOLD_BEFORE_OUTLIERS}. "
                            f"Will add interactive points only for the outlier values."
                        )

                        samples = list(value_by_sample.keys())
                        values = list(value_by_sample.values())
                        outlier_statuses = find_outliers(
                            values,
                            minval=header.get("dmin"),
                            maxval=header.get("dmax"),
                            metric=header["title"],
                        )
                        logger.debug(
                            f"Violin for '{header['title']}': found {np.count_nonzero(outlier_statuses)} outliers"
                        )
                        ds.outliers_by_metric[metric] = [
                            samples[idx] for idx in range(len(samples)) if outlier_statuses[idx]
                        ]
                        ds.show_only_outliers = True

                if values_are_numeric and not values_are_integer and "tt_decimals" in header:
                    header["hoverformat"] = f".{header['tt_decimals']}f"
                    del header["tt_decimals"]
                if "modify" in header and callable(header["modify"]):
                    del header["modify"]  # To keep the data JSON-serializable

                ds.values_by_sample_by_metric[metric] = value_by_sample
                all_samples.update(set(list(value_by_sample.keys())))

            ds.all_samples = sorted(all_samples)
            return ds

    def __init__(
        self,
        list_of_values_by_sample_by_metric: List[Dict[str, Dict[str, Union[List[int], List[float], List[str]]]]],
        list_of_header_by_metric: List[Dict[str, Dict]],
        pconfig: Dict,
        dt: Optional[DataTable] = None,
        show_table_by_default: bool = False,
    ):
        assert len(list_of_values_by_sample_by_metric) == len(list_of_header_by_metric)
        assert len(list_of_values_by_sample_by_metric) > 0

        super().__init__(PlotType.VIOLIN, pconfig, len(list_of_values_by_sample_by_metric))

        self.dt = dt
        self.show_table_by_default = dt is not None and show_table_by_default
        self.show_table = dt is not None

        self.datasets: List[ViolinPlot.Dataset] = [
            ViolinPlot.Dataset.create(ds, values_by_sample_by_metric, headers_by_metric)
            for ds, values_by_sample_by_metric, headers_by_metric in zip(
                self.datasets,
                list_of_values_by_sample_by_metric,
                list_of_header_by_metric,
            )
        ]

        self.show_only_outliers = any(ds.show_only_outliers for ds in self.datasets)

        # If the number of samples is high:
        # - do not add a table
        # - plot a Violin in Python, and serialise the figure instead of the datasets
        self.n_samples = max(len(ds.all_samples) for ds in self.datasets)
        self.serialize_figure = False
        if self.n_samples >= config.max_table_rows:
            logger.debug(
                f"The number of samples exceeds the threshold: {self.n_samples} > {config.max_table_rows}. "
                "Will build a violin plot instead of a general stats table"
            )
            self.show_table = False

        self.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            fillcolor="grey",
            line={"width": 2, "color": "grey"},
            opacity=0.5,
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

        self.layout.update(
            margin=dict(
                pad=0,
                # t=43,
                b=40,
            ),
            violingap=0,
            xaxis=dict(
                tickfont=dict(size=9, color="rgba(0,0,0,0.5)"),
            ),
            yaxis=dict(
                automargin=True,
            ),
            grid=dict(
                columns=1,
                roworder="top to bottom",
                ygap=0.4,
            ),
        )

    @staticmethod
    def from_dt(dt: DataTable, show_table_by_default=False) -> "ViolinPlot":
        values_by_sample_by_metric = dict()
        header_by_metric = dict()

        for idx, k, header in dt.get_headers_in_order():
            rid = header["rid"]
            header_by_metric[rid] = {
                "namespace": header["namespace"],
                "title": header["title"],
                "description": header["description"],
                "dmax": header.get("dmax"),
                "dmin": header.get("dmin"),
                "suffix": header.get("suffix", ""),
                "color": header.get("colour", header.get("color")),
                "hidden": header.get("hidden"),
                "modify": header.get("modify"),
                "tt_decimals": header.get("decimalPlaces", 2),
            }
            values_by_sample_by_metric[rid] = dict()
            for s_name, val_by_metric in dt.data[idx].items():
                if k in val_by_metric:
                    values_by_sample_by_metric[rid][s_name] = val_by_metric[k]

        # If all colors are the same, remove them
        if len(set([v["color"] for v in header_by_metric.values()])) == 1:
            for v in header_by_metric.values():
                v.pop("color", None)

        return ViolinPlot(
            [values_by_sample_by_metric],
            [header_by_metric],
            pconfig=dt.pconfig,
            dt=dt,
            show_table_by_default=show_table_by_default,
        )

    @staticmethod
    def tt_label() -> str:
        return ": %{x}"

    def add_to_report(self, report) -> str:
        violin_html = super().add_to_report(report)

        warning = ""
        if self.show_table_by_default and not self.show_table:
            warning = (
                f'<p class="text-muted" id="table-violin-info-{self.id}">'
                + '<span class="glyphicon glyphicon-exclamation-sign" '
                + 'title="An interactive table is not available because of the large number of samples. '
                + "A violin plot is generated instead, showing density of values for each metric, as "
                + "well as hoverable points for outlier samples in each metric. "
                + 'See http://multiqc.info/docs/#tables--beeswarm-plots"'
                + f' data-toggle="tooltip"></span> Showing {self.n_samples} samples.</p>'
            )
        elif not self.show_table:
            warning = (
                f'<p class="text-muted" id="table-violin-info-{self.id}">'
                + '<span class="glyphicon glyphicon-exclamation-sign" '
                + 'title="An interactive table is not available because of the large number of samples. '
                + "The violin plot displays hoverable points only for outlier samples in each metric, "
                + "and the hiding/highlighting functionality through the toolbox only works for outliers "
                + 'See http://multiqc.info/docs/#tables--beeswarm-plots"'
                + f' data-toggle="tooltip"></span> Showing {self.n_samples} samples.</p>'
            )

        html = (
            f"<div id='mqc-violin-{self.id}' style='{'display: none;' if (self.show_table and self.show_table_by_default) else ''}'>"
            + f"{warning}"
            + f"{violin_html}"
            + "</div>"
        )
        if self.dt:
            table_html, configuration_modal = make_table(self.dt, violin_switch=True)
            if self.show_table:
                html += f"<div id='mqc-table-{self.id}' style='{'display: none;' if not self.show_table_by_default else ''}'>{table_html}</div>"
            html += configuration_modal
        return html

    def dump_for_javascript(self):
        """Serialise the data to pick up in Plotly-JS"""
        d = super().dump_for_javascript()
        d.update(
            {
                "scatter_trace_params": self.scatter_trace_params,
                "show_only_outliers": any(ds.show_only_outliers for ds in self.datasets),
                "static": self.show_table and self.show_table_by_default,
            }
        )
        return d

    def buttons(self) -> []:
        """Add a control panel to the plot"""
        buttons = []
        if not self.flat and any(len(ds.metrics) > 1 for ds in self.datasets):
            buttons.append(
                self._btn(
                    cls="mqc_table_configModal_btn",
                    label="<span class='glyphicon glyphicon-th'></span> Configure columns",
                    data_attrs={"toggle": "modal", "target": f"#{self.id}_configModal"},
                )
            )
        if self.show_table:
            buttons.append(
                self._btn(
                    cls="mqc-violin-to-table",
                    label="<span class='glyphicon glyphicon-th-list'></span> Table",
                )
            )

        return buttons + super().buttons()

    def create_figure(
        self,
        layout: go.Layout,
        dataset: Dataset,
        is_log=False,
        is_pct=False,
        add_scatter=True,
    ):
        """
        Create a Plotly figure for a dataset
        """
        metrics = [m for m in dataset.metrics if not dataset.header_by_metric[m].get("hidden", False)]

        layout = copy.deepcopy(layout)
        layout.height = 70 * len(metrics) + 50
        layout.grid.rows = len(metrics)
        layout.grid.subplots = [[(f"x{i + 1}y{i + 1}" if i > 0 else "xy")] for i in range(len(metrics))]

        for metric_idx, metric in enumerate(metrics):
            header = dataset.header_by_metric[metric]
            layout[f"xaxis{metric_idx + 1}"] = {
                "gridcolor": layout["xaxis"]["gridcolor"],
                "zerolinecolor": layout["xaxis"]["zerolinecolor"],
                "tickfont": copy.deepcopy(layout["xaxis"]["tickfont"]),
            }
            layout[f"xaxis{metric_idx + 1}"].update(header.get("xaxis", {}))
            layout[f"yaxis{metric_idx + 1}"] = {
                "gridcolor": layout["yaxis"]["gridcolor"],
                "zerolinecolor": layout["yaxis"]["zerolinecolor"],
                "automargin": True,
            }
            title = header["title"] + "  "
            if header.get("namespace"):
                title = f"{header['namespace']}  <br>" + title
            layout[f"yaxis{metric_idx + 1}"].update(
                {
                    "tickmode": "array",
                    "tickvals": [metric_idx],
                    "ticktext": [title],
                }
            )
            if header.get("color"):
                layout[f"yaxis{metric_idx + 1}"]["tickfont"] = {
                    "color": f"rgb({header['color']})",
                }

        layout.xaxis = layout["xaxis1"]
        layout.yaxis = layout["yaxis1"]

        layout.showlegend = False
        fig = go.Figure(layout=layout)

        for metric_idx, metric in enumerate(metrics):
            header = dataset.header_by_metric[metric]
            values_by_sample = dataset.values_by_sample_by_metric[metric]
            outliers = dataset.outliers_by_metric.get(metric, [])

            params = copy.deepcopy(self.trace_params)
            color = header.get("color")
            if color:
                params["fillcolor"] = f"rgba({color},1)"
                params["line"]["color"] = f"rgba({color},1)"

            axis_key = "" if metric_idx == 0 else str(metric_idx + 1)
            fig.add_trace(
                go.Violin(
                    x=list(values_by_sample.values()),
                    name=metric_idx,
                    text=list(values_by_sample.keys()),
                    xaxis=f"x{axis_key}",
                    yaxis=f"y{axis_key}",
                    **params,
                ),
            )
            if add_scatter:
                for sample, value in values_by_sample.items():
                    scatter_params = copy.deepcopy(self.scatter_trace_params)
                    scatter_params["showlegend"] = False
                    if dataset.show_only_outliers and sample not in outliers:
                        continue
                    fig.add_trace(
                        go.Scatter(
                            x=[value],
                            y=[metric_idx],
                            text=[sample],
                            xaxis=f"x{axis_key}",
                            yaxis=f"y{axis_key}",
                            **scatter_params,
                        ),
                    )
        return fig

    def save_data_file(self, data: BaseDataset) -> None:
        pass


def find_outliers(
    values: Union[List[int], List[float]],
    top_n: Optional[int] = None,
    z_cutoff: float = 2.0,
    minval: Optional[Union[float, int]] = None,
    maxval: Optional[Union[float, int]] = None,
    metric: Optional[str] = None,
) -> np.array:
    """
    If `n` is defined, find `n` most outlying points in a list.
    Otherwise, find outliers with a Z-score above `z_cutoff`.

    Return a list of booleans, indicating if the value with this index is an outlier.

    `minval` and/or `maxval` can be added to include them into the outlier detection array.
    This is a trick to avoid the following problem: for example, percentage values are
    clustered around 0, but differ in e.g. 3+rd decimal place. In this case, the outliers
    will be nearly no different from the other values. If we add "100%" as an artificial
    additional value, none of those near-zero clustered values won't be called outliers.
    """
    if len(values) == 0 or (top_n is not None and top_n <= 0):
        return np.zeros(len(values), dtype=bool)

    added_values = []
    if minval is not None:
        added_values.append(minval)
    if maxval is not None:
        added_values.append(maxval)
    values = np.array(values + added_values)

    # Calculate the mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    if std_dev == 0:
        logger.warning(f"All {len(values)} points have the same values" + (f", metric: '{metric}'" if metric else ""))
        return np.zeros(len(values), dtype=bool)

    # Calculate Z-scores (measures of "outlyingness")
    z_scores = np.abs((values - mean) / std_dev)

    # Get indices of the top N outliers
    outlier_status = np.zeros(len(values), dtype=bool)
    if top_n:
        outlier_status[np.argsort(z_scores)[-top_n:]] = True
    else:
        indices = np.where(z_scores > z_cutoff)[0]
        while len(indices) <= len(added_values) and z_cutoff > 1.0:
            new_z_cutoff = z_cutoff - 0.2
            logger.warning(
                f"No outliers found with Z-score cutoff {z_cutoff:.1f}, trying a lower cutoff: {new_z_cutoff:.1f}"
                + (f", metric: '{metric}'" if metric else "")
            )
            z_cutoff = new_z_cutoff
            indices = np.where(z_scores > z_cutoff)[0]
        outlier_status[indices] = True

    if added_values:
        outlier_status = outlier_status[: -len(added_values)]
    return outlier_status
