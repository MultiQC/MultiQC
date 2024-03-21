import dataclasses
import logging
from typing import Dict, List, Union, Any, Optional
import copy

import math
import numpy as np
import plotly.graph_objects as go

from multiqc.utils import config, util_functions
from multiqc.plots.table_object import DataTable
from multiqc.plots.plotly.plot import Plot, PlotType, BaseDataset
from multiqc.plots.plotly.table import make_table

logger = logging.getLogger(__name__)


def plot(dts: List[DataTable], show_table_by_default=False) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    """
    p = ViolinPlot.from_dt(dts, show_table_by_default)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class Dataset(BaseDataset):
    metrics: List[str]
    header_by_metric: Dict[str, Dict[str, Any]]
    violin_values_by_sample_by_metric: Dict[str, Dict[str, Union[List[int], List[float], List[str]]]]
    scatter_values_by_sample_by_metric: Dict[str, Dict[str, Union[List[int], List[float], List[str]]]]
    all_samples: List[str]  # unique list of all samples in this dataset
    scatter_trace_params: Dict[str, Any]

    @staticmethod
    def create(
        dataset: BaseDataset,
        values_by_sample_by_metric: Dict[str, Dict[str, Union[List[int], List[float], List[str]]]],
        header_by_metric: Dict[str, Dict[str, Any]],
    ) -> "Dataset":
        ds = Dataset(
            **dataset.__dict__,
            metrics=[],
            header_by_metric=header_by_metric,
            violin_values_by_sample_by_metric=dict(),
            scatter_values_by_sample_by_metric=dict(),
            all_samples=[],
            scatter_trace_params=dict(),
        )

        all_samples = set()
        for metric, header in header_by_metric.items():
            # Add Plotly-specific parameters to the header
            xaxis = {"ticksuffix": header.get("suffix")}
            header["xaxis"] = xaxis
            header["show_points"] = True
            header["show_only_outliers"] = False

            # Take non-empty values for the violin
            value_by_sample = values_by_sample_by_metric[metric]
            value_by_sample = {s: v for s, v in value_by_sample.items() if v is not None and str(v).strip() != ""}

            # Parse values to numbers if possible
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
                    # noinspection PyBroadException
                    try:
                        value_by_sample[s] = header["modify"](v)
                    except Exception:  # User-provided modify function can raise any exception
                        pass

            if not value_by_sample:
                logger.debug(f"No non-empty values found for metric: {header['title']}")
                if "modify" in header and callable(header["modify"]):
                    del header["modify"]  # To keep the data JSON-serializable
                continue

            values_are_numeric = all(isinstance(v, (int, float)) for v in value_by_sample.values())
            values_are_integer = all(isinstance(v, int) for v in value_by_sample.values())
            if values_are_numeric:
                # Remove NaN and Inf values
                value_by_sample = {s: v for s, v in value_by_sample.items() if np.isfinite(v)}
                if not value_by_sample:
                    logger.warning(f"All values are NaN or Inf for metric: {header['title']}")
                    if "modify" in header and callable(header["modify"]):
                        del header["modify"]  # To keep the data JSON-serializable
                    continue

            header["show_points"] = len(value_by_sample) <= config.violin_min_threshold_no_points
            header["show_only_outliers"] = len(value_by_sample) > config.violin_min_threshold_outliers

            if values_are_numeric:
                # Calculate range
                xmin = header.get("dmin")
                xmax = header.get("dmax")
                tickvals = None
                if xmin == xmax == 0:  # Plotly will modify the 0:0 range to -1:1, and we want to keep it 0-centered
                    xmax = 1
                    tickvals = [0]
                xaxis["tickvals"] = tickvals
                if xmin is not None and xmax is not None:
                    if header["show_points"]:
                        # add extra padding to avoid clipping the points at range limits
                        xmin -= (xmax - xmin) * 0.005
                        xmax += (xmax - xmin) * 0.005
                    xaxis["range"] = [xmin, xmax]

            if not header["show_points"]:  # Do not add any interactive points
                scatter_values_by_sample = {}
            elif not header["show_only_outliers"]:
                scatter_values_by_sample = {}  # will use the violin values
            else:
                if not values_are_numeric:
                    logger.debug(
                        f"Violin for '{header['title']}': sample number is {len(value_by_sample)} > {header['show_only_outliers']}. "
                        f"As values are not numeric, will not add any interactive points."
                    )
                    scatter_values_by_sample = {}
                else:
                    logger.debug(
                        f"Violin for '{header['title']}': sample number is {len(value_by_sample)} > {header['show_only_outliers']}. "
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
                    logger.debug(f"Violin for '{header['title']}': found {np.count_nonzero(outlier_statuses)} outliers")
                    scatter_values_by_sample = {
                        samples[idx]: values[idx] for idx in range(len(samples)) if outlier_statuses[idx]
                    }

            ds.scatter_values_by_sample_by_metric[metric] = scatter_values_by_sample

            # Now sort and downsample values to keep max 2000 points for each metric
            violin_values_by_sample = value_by_sample
            max_violin_points = config.violin_downsample_after
            if max_violin_points is not None and len(violin_values_by_sample) > max_violin_points:
                logger.debug(
                    f"Violin for '{header['title']}': sample number is {len(violin_values_by_sample)}. "
                    f"Will downsample to max {max_violin_points} points."
                )
                samples = list(violin_values_by_sample.keys())
                values = list(violin_values_by_sample.values())
                indices = np.argsort(values)
                indices = indices[:: int(math.ceil(len(indices) / max_violin_points))]
                violin_values_by_sample = {samples[idx]: values[idx] for idx in indices}

            ds.violin_values_by_sample_by_metric[metric] = violin_values_by_sample

            # Clean up the header
            if values_are_numeric and not values_are_integer and "tt_decimals" in header:
                header["hoverformat"] = f".{header['tt_decimals']}f"
                del header["tt_decimals"]
            if "modify" in header and callable(header["modify"]):
                del header["modify"]  # To keep the data JSON-serializable

            all_samples.update(set(list(scatter_values_by_sample.keys())))
            all_samples.update(set(list(violin_values_by_sample.keys())))
            ds.metrics.append(metric)

        ds.all_samples = sorted(all_samples)

        ds.trace_params.update(
            orientation="h",
            box={"visible": True},
            meanline={"visible": True},
            fillcolor="#b5b5b5",
            line={"width": 2, "color": "#b5b5b5"},
            opacity=0.5,
            points=False,  # Don't show points, we'll add them manually
            # The hover information is useful, but the formatting is ugly and not
            # configurable as far as I can see. Also, it's not possible to disable it,
            # so setting it to "points" as we don't show points, so it's effectively
            # disabling it.
            hoveron="points",
        )

        # If all violins are grey, make the dots blue to make it more clear that it's interactive
        # if some violins are color-coded, make the dots black to make them less distracting
        marker_color = "black" if any(h.get("color") for h in ds.header_by_metric.values()) else "#0b79e6"
        ds.scatter_trace_params = {
            "mode": "markers",
            "marker": {
                "size": 4,
                "color": marker_color,
            },
            "showlegend": False,
            "hovertemplate": ds.trace_params["hovertemplate"],
            "hoverlabel": {"bgcolor": "white"},
        }
        return ds

    def create_figure(
        self,
        layout: go.Layout,
        is_log=False,
        is_pct=False,
        add_scatter=True,
    ):
        """
        Create a Plotly figure for a dataset
        """
        metrics = [m for m in self.metrics if not self.header_by_metric[m].get("hidden", False)]

        layout = copy.deepcopy(layout)
        layout.grid.rows = len(metrics)
        layout.grid.subplots = [[(f"x{i + 1}y{i + 1}" if i > 0 else "xy")] for i in range(len(metrics))]
        layout.height = ViolinPlot.VIOLIN_HEIGHT * len(metrics) + ViolinPlot.EXTRA_HEIGHT

        for metric_idx, metric in enumerate(metrics):
            header = self.header_by_metric[metric]

            layout[f"xaxis{metric_idx + 1}"] = {
                "automargin": layout["xaxis"]["automargin"],
                "color": layout["xaxis"]["color"],
                "gridcolor": layout["xaxis"]["gridcolor"],
                "zerolinecolor": layout["xaxis"]["zerolinecolor"],
                "hoverformat": layout["xaxis"]["hoverformat"],
                "tickfont": copy.deepcopy(layout["xaxis"]["tickfont"]),
            }
            layout[f"xaxis{metric_idx + 1}"].update(header.get("xaxis", {}))
            layout[f"yaxis{metric_idx + 1}"] = {
                "automargin": layout["yaxis"]["automargin"],
                "color": layout["yaxis"]["color"],
                "gridcolor": layout["yaxis"]["gridcolor"],
                "zerolinecolor": layout["yaxis"]["zerolinecolor"],
                "hoverformat": layout["yaxis"]["hoverformat"],
                "tickfont": copy.deepcopy(layout["yaxis"]["tickfont"]),
            }

            padding = "  "  # otherwise the labels will stick too close to the axis
            title = header["title"] + padding
            if header.get("namespace"):
                title = f"{header['namespace']}{padding}<br>" + title
            layout[f"yaxis{metric_idx + 1}"].update(
                {
                    "tickmode": "array",
                    "tickvals": [metric_idx],
                    "ticktext": [title],
                }
            )

            if "hoverformat" in header:
                layout[f"xaxis{metric_idx + 1}"]["hoverformat"] = header["hoverformat"]

            if header.get("color"):
                layout[f"yaxis{metric_idx + 1}"]["tickfont"] = {
                    "color": f"rgb({header['color']})",
                }

        layout["xaxis"] = layout["xaxis1"]
        layout["yaxis"] = layout["yaxis1"]

        fig = go.Figure(layout=layout)

        violin_values_by_sample_by_metric = self.violin_values_by_sample_by_metric

        for metric_idx, metric in enumerate(metrics):
            header = self.header_by_metric[metric]
            params = copy.deepcopy(self.trace_params)
            color = header.get("color")
            if color:
                params["fillcolor"] = f"rgb({color})"
                params["line"]["color"] = f"rgb({color})"

            violin_values_by_sample = violin_values_by_sample_by_metric[metric]
            axis_key = "" if metric_idx == 0 else str(metric_idx + 1)
            fig.add_trace(
                go.Violin(
                    x=list(violin_values_by_sample.values()),
                    name=metric_idx,
                    text=list(violin_values_by_sample.keys()),
                    xaxis=f"x{axis_key}",
                    yaxis=f"y{axis_key}",
                    **params,
                ),
            )

            if add_scatter and header["show_points"]:
                if header["show_only_outliers"]:
                    scatter_values_by_sample = self.scatter_values_by_sample_by_metric[metric]
                else:
                    scatter_values_by_sample = violin_values_by_sample
                scatter_params = copy.deepcopy(self.scatter_trace_params)
                for sample, value in scatter_values_by_sample.items():
                    # add vertical jitter (not working in python version currently)
                    y = float(metric_idx)
                    # y += random.uniform(-0.2, 0.2)
                    # y += random.random() * 0.3 - 0.3 / 2
                    fig.add_trace(
                        go.Scatter(
                            x=[value],
                            y=[y],
                            text=[sample],
                            xaxis=f"x{axis_key}",
                            yaxis=f"y{axis_key}",
                            **scatter_params,
                        ),
                    )
        return fig


class ViolinPlot(Plot):
    VIOLIN_HEIGHT = 70  # single violin height
    EXTRA_HEIGHT = 63  # extra space for the title and footer

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

        super().__init__(
            PlotType.VIOLIN,
            pconfig,
            n_datasets=len(list_of_values_by_sample_by_metric),
            id=dt.id if dt else None,
        )

        self.dt = dt
        self.no_violin = pconfig.get("no_violin", pconfig.get("no_beeswarm", False))
        self.show_table = dt is not None
        self.show_table_by_default = show_table_by_default or self.no_violin
        if dt:  # to make it different from the violin id
            dt.id = "table-" + dt.id

        self.datasets: List[Dataset] = [
            Dataset.create(ds, values_by_sample_by_metric, headers_by_metric)
            for ds, values_by_sample_by_metric, headers_by_metric in zip(
                self.datasets,
                list_of_values_by_sample_by_metric,
                list_of_header_by_metric,
            )
        ]

        # Violin-specific layout parameters
        self.layout.update(
            margin=dict(
                pad=0,
                b=40,
            ),
            xaxis=dict(
                # so Plotly doesn't try to fit the ticks the on the most bottom violin,
                # and squish the other violins
                automargin=False,
                tickfont=dict(size=9, color="rgba(0,0,0,0.5)"),
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
            ),
            yaxis=dict(
                tickfont=dict(size=9, color="rgba(0,0,0,0.5)"),
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
            ),
            violingap=0,
            grid=dict(
                columns=1,
                roworder="top to bottom",
                ygap=0.4,
            ),
            showlegend=False,
        )

        # If the number of samples is high:
        # - do not add a table
        # - plot a Violin in Python, and serialise the figure instead of the datasets
        self.n_samples = max(len(ds.all_samples) for ds in self.datasets)
        self.serialize_figure = False
        if self.n_samples > config.max_table_rows and not self.no_violin:
            self.show_table = False
            if self.show_table_by_default:
                logger.debug(
                    f"Table '{self.id}': sample number {self.n_samples} > {config.max_table_rows}, "
                    "Will render only a violin plot instead of the table"
                )

    @staticmethod
    def from_dt(dts: List[DataTable], show_table_by_default=False) -> "ViolinPlot":
        list_values_by_sample_by_metric = []
        list_header_by_metric = []

        for dt in dts:
            list_values_by_sample_by_metric.append(dict())
            list_header_by_metric.append(dict())

            for idx, k, header in dt.get_headers_in_order():
                rid = header["rid"]
                list_header_by_metric[-1][rid] = {
                    "namespace": header["namespace"],
                    "title": header["title"],
                    "description": header["description"],
                    "dmax": header.get("dmax"),
                    "dmin": header.get("dmin"),
                    "suffix": header.get("suffix", ""),
                    "color": header.get("colour", header.get("color")),
                    "hidden": header.get("hidden"),
                    "modify": header.get("modify"),
                    "tt_decimals": header.get("tt_decimals", header.get("decimalPlaces", 2)),
                }
                list_values_by_sample_by_metric[-1][rid] = dict()
                for s_name, val_by_metric in dt.data[idx].items():
                    if k in val_by_metric:
                        list_values_by_sample_by_metric[-1][rid][s_name] = val_by_metric[k]

            # If all colors are the same, remove them
            if len(set([v["color"] for v in list_header_by_metric[-1].values()])) == 1:
                for v in list_header_by_metric[-1].values():
                    v.pop("color", None)

            # If all namespaces are the same as well, remove them too (usually they follow the colors pattern)
            if len(set([v["namespace"] for v in list_header_by_metric[-1].values()])) == 1:
                for v in list_header_by_metric[-1].values():
                    v.pop("namespace", None)

        return ViolinPlot(
            list_values_by_sample_by_metric,
            list_header_by_metric,
            pconfig=dts[0].pconfig,
            dt=dts[0],
            show_table_by_default=show_table_by_default,
        )

    @staticmethod
    def tt_label() -> str:
        return ": %{x}"

    def add_to_report(self, report) -> str:
        warning = ""
        if self.show_table_by_default and not self.show_table:
            warning = (
                f'<p class="text-muted" id="table-violin-info-{self.id}">'
                + '<span class="glyphicon glyphicon-exclamation-sign" '
                + 'title="An interactive table is not available because of the large number of samples. '
                + "A violin plot is generated instead, showing density of values for each metric, as "
                + 'well as hoverable points for outlier samples in each metric."'
                + f' data-toggle="tooltip"></span> Showing {self.n_samples} samples.</p>'
            )
        elif not self.show_table:
            warning = (
                f'<p class="text-muted" id="table-violin-info-{self.id}">'
                + '<span class="glyphicon glyphicon-exclamation-sign" '
                + 'title="An interactive table is not available because of the large number of samples. '
                + "The violin plot displays hoverable points only for outlier samples in each metric, "
                + 'and the hiding/highlighting functionality through the toolbox only works for outliers"'
                + f' data-toggle="tooltip"></span> Showing {self.n_samples} samples.</p>'
            )

        if not self.show_table:
            # Show violin alone.
            # Note that "no_violin" will be ignored here as we need to render _something_. The only case it can
            # happen if violin.plot() is called directly, and "no_violin" is passed, which doesn't make sense.
            html = warning + super().add_to_report(report)
        elif self.no_violin:
            # Show table alone
            table_html, configuration_modal = make_table(self.dt)
            html = warning + table_html + configuration_modal
        else:
            # Render both, add a switch between table and violin
            table_html, configuration_modal = make_table(self.dt, violin_id=self.id)
            violin_html = super().add_to_report(report)

            violin_visibility = "style='display: none;'" if self.show_table_by_default else ""
            html = f"<div id='mqc_violintable_wrapper_{self.id}' {violin_visibility}>{warning}{violin_html}</div>"

            table_visibility = "style='display: none;'" if not self.show_table_by_default else ""
            html += f"<div id='mqc_violintable_wrapper_{self.dt.id}' {table_visibility}>{table_html}</div>"

            html += configuration_modal

        return html

    def dump_for_javascript(self):
        """Serialise the data to pick up in Plotly-JS"""
        d = super().dump_for_javascript()
        d.update(
            {
                "static": self.show_table and self.show_table_by_default,
                "violin_height": ViolinPlot.VIOLIN_HEIGHT,
                "extra_height": ViolinPlot.EXTRA_HEIGHT,
            }
        )
        return d

    def buttons(self) -> []:
        """Add a control panel to the plot"""
        buttons = []
        if not self.flat and any(len(ds.metrics) > 1 for ds in self.datasets) and self.dt is not None:
            buttons.append(
                self._btn(
                    cls="mqc_table_configModal_btn",
                    label="<span class='glyphicon glyphicon-th'></span> Configure columns",
                    data_attrs={"toggle": "modal", "target": f"#{self.dt.id}_configModal"},
                )
            )
        if self.show_table:
            buttons.append(
                self._btn(
                    cls="mqc-violin-to-table",
                    label="<span class='glyphicon glyphicon-th-list'></span> Table",
                    data_attrs={"table-id": self.dt.id, "violin-id": self.id},
                )
            )

        return buttons + super().buttons()

    def save_data_file(self, dataset: Dataset) -> None:
        data = {}
        for metric in dataset.metrics:
            values_by_sample = dataset.violin_values_by_sample_by_metric[metric]
            title = dataset.header_by_metric[metric]["title"]
            for sample, value in values_by_sample.items():
                data.setdefault(sample, {})[title] = value

        util_functions.write_data_file(data, dataset.uid)


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
        logger.debug(f"All {len(values)} points have the same values" + (f", metric: '{metric}'" if metric else ""))
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
            logger.debug(
                f"No outliers found with Z-score cutoff {z_cutoff:.1f}, trying a lower cutoff: {new_z_cutoff:.1f}"
                + (f", metric: '{metric}'" if metric else "")
            )
            z_cutoff = new_z_cutoff
            indices = np.where(z_scores > z_cutoff)[0]
        outlier_status[indices] = True

    if added_values:
        outlier_status = outlier_status[: -len(added_values)]
    return outlier_status
