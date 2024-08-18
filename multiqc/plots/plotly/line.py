import io
import logging
import os
import random
from typing import Dict, List, Union, Tuple, Optional, Literal, Any, Mapping, TypeVar, Generic

import math
import plotly.graph_objects as go  # type: ignore
from plotly.graph_objs.layout.shape import Label  # type: ignore
from pydantic import Field, BaseModel

from multiqc.plots.plotly.plot import PlotType, BaseDataset, Plot, PConfig
from multiqc import config, report
from multiqc.utils.util_functions import update_dict
from multiqc.validation import ValidatedConfig, add_validation_warning

logger = logging.getLogger(__name__)


KeyT = TypeVar("KeyT", int, str, float)
ValueT = TypeVar("ValueT", int, str, float, None)
XToYDictT = Mapping[KeyT, ValueT]
DatasetT = Mapping[str, XToYDictT]


class Marker(BaseModel):
    symbol: Optional[str] = None
    color: Optional[str] = None
    line_color: Optional[str] = None
    fill_color: Optional[str] = None
    width: int = 1


class Series(ValidatedConfig, Generic[KeyT, ValueT]):
    name: str = Field(default_factory=lambda: f"series-{random.randint(1000000, 9999999)}")
    pairs: List[Tuple[KeyT, ValueT]]
    color: Optional[str] = None
    width: int = 2
    dash: Optional[str] = None
    showlegend: bool = True
    marker: Optional[Marker] = None

    def __init__(self, **data):
        if "dashStyle" in data:
            add_validation_warning(
                [LinePlotConfig, Series], "'dashStyle' field is deprecated. Please use 'dash' instead"
            )
            data["dash"] = convert_dash_style(Series, data.pop("dashStyle"))
        elif "dash" in data:
            data["dash"] = convert_dash_style(Series, data["dash"])

        tuples: List[Tuple[KeyT, ValueT]] = []
        if "data" in data:
            add_validation_warning([LinePlotConfig, Series], "'data' field is deprecated. Please use 'pairs' instead")
        for p in data.pop("data") if "data" in data else data.get("pairs", []):
            if isinstance(p, list):
                tuples.append(tuple(p))
            else:
                tuples.append(p)
        data["pairs"] = tuples

        super().__init__(**data, _parent_class=LinePlotConfig)


class FlatLine(ValidatedConfig):
    """
    Extra X=const or Y=const line added to the plot
    """

    value: Union[float, int]
    colour: Optional[str] = Field(None, deprecated="color")
    color: Optional[str] = None
    width: int = 2
    dashStyle: Optional[str] = Field(None, deprecated="dash")
    dash: Optional[str] = None
    label: Optional[Union[str, Dict]] = None

    @classmethod
    def parse_label(cls, value):
        if isinstance(value, dict):
            add_validation_warning(
                [LinePlotConfig, FlatLine],
                "Line plot's x_lines or y_lines 'label' field is expected to be a string. "
                "Other fields other than 'text' are deprecated and will be ignored",
            )
            return value["text"]
        return value

    def __init__(self, **data):
        if "dashStyle" in data:
            data["dash"] = convert_dash_style(FlatLine, data.pop("dashStyle"))
        if "dash" in data:
            data["dash"] = convert_dash_style(FlatLine, data["dash"])
        super().__init__(**data, _parent_class=LinePlotConfig)


class LineBand(ValidatedConfig):
    """
    Extra X1-X2 or Y1-Y2 band added to the plot
    """

    from_: Union[float, int]
    to: Union[float, int]
    color: Optional[str] = None

    def __init__(self, **data):
        super().__init__(**data, _parent_class=LinePlotConfig)


SeriesConf = Union[Series, Dict]


class LinePlotConfig(PConfig):
    xlab: Optional[str] = None
    ylab: Optional[str] = None
    categories: bool = False
    smooth_points: Optional[int] = 500
    smooth_points_sumcounts: Union[bool, List[bool], None] = None
    extra_series: Optional[Union[SeriesConf, List[SeriesConf], List[List[SeriesConf]]]] = None
    xMinRange: Optional[Union[float, int]] = Field(None, deprecated="x_minrange")
    yMinRange: Optional[Union[float, int]] = Field(None, deprecated="y_minrange")
    x_minrange: Optional[Union[float, int]] = None
    y_minrange: Optional[Union[float, int]] = None
    xPlotBands: Optional[List[LineBand]] = Field(None, deprecated="x_bands")
    yPlotBands: Optional[List[LineBand]] = Field(None, deprecated="y_bands")
    xPlotLines: Optional[List[FlatLine]] = Field(None, deprecated="x_lines")
    yPlotLines: Optional[List[FlatLine]] = Field(None, deprecated="y_lines")
    x_bands: Optional[List[LineBand]] = None
    y_bands: Optional[List[LineBand]] = None
    x_lines: Optional[List[FlatLine]] = None
    y_lines: Optional[List[FlatLine]] = None
    style: Optional[Literal["lines", "lines+markers"]] = None
    hide_zero_cats: Optional[bool] = Field(False, deprecated="hide_empty")
    hide_empty: bool = False
    colors: Dict[str, str] = {}

    @classmethod
    def parse_extra_series(cls, data: Union[SeriesConf, List[SeriesConf], List[List[SeriesConf]]]):
        if isinstance(data, list):
            if isinstance(data[0], list):
                return [[Series(**d) if isinstance(d, dict) else d for d in ds] for ds in data]
            return [Series(**d) if isinstance(d, dict) else d for d in data]
        return Series(**data) if isinstance(data, dict) else data

    @classmethod
    def parse_x_bands(cls, data):
        return [LineBand(**d) for d in data]

    @classmethod
    def parse_y_bands(cls, data):
        return [LineBand(**d) for d in data]

    @classmethod
    def parse_x_lines(cls, data):
        return [FlatLine(**d) for d in data]

    @classmethod
    def parse_y_lines(cls, data):
        return [FlatLine(**d) for d in data]


def plot(lists_of_lines: List[List[Series]], pconfig: LinePlotConfig) -> "LinePlot":
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param lists_of_lines: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """

    # if self.n_samples >= config.max_table_rows:
    # Get a line of median values at each point with an interval of max and
    # min values
    # Create a violin of median values in each sample, showing dots for outliers
    # Clicking on a dot of a violin will show the line plot for that sample

    return create(pconfig, lists_of_lines)


class Dataset(BaseDataset):
    lines: List[Series]

    @staticmethod
    def create(
        dataset: BaseDataset,
        lines: List[Series],
        pconfig: LinePlotConfig,
    ) -> "Dataset":
        dataset = Dataset(**dataset.model_dump(), lines=lines)

        # Prevent Plotly-JS from parsing strings as numbers
        if pconfig.categories or dataset.dconfig.get("categories"):
            dataset.layout["xaxis"]["type"] = "category"

        if pconfig.style is not None:
            mode = pconfig.style
        else:
            num_data_points = sum(len(x.pairs) for x in lines)
            if num_data_points < config.lineplot_number_of_points_to_hide_markers:
                mode = "lines+markers"
            else:
                mode = "lines"

        dataset.trace_params.update(
            mode=mode,
            line={"width": 2},
        )
        if mode == "lines+markers":
            dataset.trace_params.update(
                line={"width": 0.6},
                marker={"size": 5},
            )
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log=False,
        is_pct=False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        if layout.showlegend is True:
            # Extra space for legend
            layout.height += len(self.lines) * 5

        fig = go.Figure(layout=layout)
        for series in self.lines:
            xs = [x[0] for x in series.pairs]
            ys = [x[1] for x in series.pairs]
            if series.dash:
                print(series)
            params: Dict[str, Any] = {
                "showlegend": series.showlegend,
                "line": {
                    "color": series.color,
                    "dash": series.dash,
                    "width": series.width,
                },
            }
            if series.marker:
                params["mode"] = "lines+markers"
                params["marker"] = {
                    "symbol": series.marker.symbol,
                    "color": series.marker.fill_color or series.marker.color or series.color,
                    "line": {
                        "width": series.marker.width,
                        "color": series.marker.line_color or series.marker.color or "black",
                    },
                }
            params = update_dict(params, self.trace_params, none_only=True)
            if len(series.pairs) == 1:
                params["mode"] = "lines+markers"  # otherwise it's invisible

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=series.name,
                    text=[series.name] * len(xs),
                    **params,
                )
            )
        return fig

    def save_data_file(self) -> None:
        y_by_x_by_sample: Dict[str, Dict] = dict()
        last_cats = None
        shared_cats = True
        for series in self.lines:
            y_by_x_by_sample[series.name] = dict()

            # Check to see if all categories are the same
            if len(series.pairs) > 0 and isinstance(series.pairs[0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in series.pairs]
                elif last_cats != [x[0] for x in series.pairs]:
                    shared_cats = False

            for i, x in enumerate(series.pairs):
                if isinstance(x, list):
                    y_by_x_by_sample[series.name][x[0]] = x[1]
                else:
                    try:
                        y_by_x_by_sample[series.name][self.dconfig["categories"][i]] = x
                    except (ValueError, KeyError, IndexError):
                        y_by_x_by_sample[series.name][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format in ["tsv", "csv"]:
            sep = "\t" if config.data_format == "tsv" else ","
            fout = ""
            for series in self.lines:
                fout += series.name + sep + "X" + sep + sep.join([str(x[0]) for x in series.pairs]) + "\n"
                fout += series.name + sep + "Y" + sep + sep.join([str(x[1]) for x in series.pairs]) + "\n"

            fn = f"{self.uid}.{config.data_format_extensions[config.data_format]}"
            fpath = os.path.join(report.data_tmp_dir(), fn)
            with io.open(fpath, "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            report.write_data_file(y_by_x_by_sample, self.uid)


class LinePlot(Plot[Dataset]):
    datasets: List[Dataset]


def create(
    pconfig: LinePlotConfig,
    lists_of_lines: List[List[Series]],
) -> "LinePlot":
    n_samples_per_dataset = [len(x) for x in lists_of_lines]

    model = Plot.initialize(
        plot_type=PlotType.LINE,
        pconfig=pconfig,
        n_samples_per_dataset=n_samples_per_dataset,
        axis_controlled_by_switches=["yaxis"],
        default_tt_label="<br>%{x}: %{y}",
    )

    # Very large legend for automatically enabled flat plot mode is not very helpful
    max_n_samples = max(len(x) for x in lists_of_lines) if len(lists_of_lines) > 0 else 0
    if pconfig.showlegend is None and max_n_samples > 250:
        model.layout.showlegend = False

    model.datasets = [Dataset.create(d, lines, pconfig) for d, lines in zip(model.datasets, lists_of_lines)]

    # Make a tooltip always show on hover over any point on plot
    model.layout.hoverdistance = -1

    y_minrange = pconfig.y_minrange
    x_minrange = pconfig.x_minrange
    y_bands = pconfig.y_bands
    x_bands = pconfig.x_bands
    x_lines = pconfig.x_lines
    y_lines = pconfig.y_lines
    if y_minrange or y_bands or y_lines:
        # We don't want the bands to affect the calculated axis range, so we
        # find the min and the max from data points, and manually set the range.
        for dataset in model.datasets:
            minval = dataset.layout["yaxis"]["autorangeoptions"]["minallowed"]
            maxval = dataset.layout["yaxis"]["autorangeoptions"]["maxallowed"]
            for line in dataset.lines:
                ys = [x[1] for x in line.pairs]
                if len(ys) > 0:
                    minval = min(ys) if minval is None else min(minval, min(ys))
                    maxval = max(ys) if maxval is None else max(maxval, max(ys))
            if maxval is not None and minval is not None:
                maxval += (maxval - minval) * 0.05
            clipmin = dataset.layout["yaxis"]["autorangeoptions"]["clipmin"]
            clipmax = dataset.layout["yaxis"]["autorangeoptions"]["clipmax"]
            if clipmin is not None and minval is not None and clipmin > minval:
                minval = clipmin
            if clipmax is not None and maxval is not None and clipmax < maxval:
                maxval = clipmax
            if y_minrange is not None and maxval is not None and minval is not None:
                maxval = max(maxval, minval + y_minrange)
            if model.layout.yaxis.type == "log":
                minval = math.log10(minval) if minval is not None and minval > 0 else None
                maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
            dataset.layout["yaxis"]["range"] = [minval, maxval]

    if not pconfig.categories and x_minrange or x_bands or x_lines:
        # same as above but for x-axis
        for dataset in model.datasets:
            minval = dataset.layout["xaxis"]["autorangeoptions"]["minallowed"]
            maxval = dataset.layout["xaxis"]["autorangeoptions"]["maxallowed"]
            for series in dataset.lines:
                xs = [x[0] for x in series.pairs]
                if len(xs) > 0:
                    minval = min(xs) if minval is None else min(minval, min(xs))
                    maxval = max(xs) if maxval is None else max(maxval, max(xs))
            clipmin = dataset.layout["xaxis"]["autorangeoptions"]["clipmin"]
            clipmax = dataset.layout["xaxis"]["autorangeoptions"]["clipmax"]
            if clipmin is not None and minval is not None and clipmin > minval:
                minval = clipmin
            if clipmax is not None and maxval is not None and clipmax < maxval:
                maxval = clipmax
            if x_minrange is not None and maxval is not None and minval is not None:
                maxval = max(maxval, minval + x_minrange)
            if model.layout.xaxis.type == "log":
                minval = math.log10(minval) if minval is not None and minval > 0 else None
                maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
            dataset.layout["xaxis"]["range"] = [minval, maxval]

    model.layout.shapes = (
        [
            dict(
                type="rect",
                y0=band.from_,
                y1=band.to,
                x0=0,
                x1=1,
                xref="paper",  # make x coords are relative to the plot paper [0,1]
                fillcolor=band.color,
                line={
                    "width": 0,
                },
                layer="below",
            )
            for band in (y_bands or [])
        ]
        + [
            dict(
                type="rect",
                x0=band.from_,
                x1=band.to,
                y0=0,
                y1=1,
                yref="paper",  # make y coords are relative to the plot paper [0,1]
                fillcolor=band.color,
                line={
                    "width": 0,
                },
                layer="below",
            )
            for band in (x_bands or [])
        ]
        + [
            dict(
                type="line",
                xref="paper",
                yref="y",
                x0=0,
                y0=line.value,
                x1=1,
                y1=line.value,
                line={
                    "width": line.width,
                    "dash": line.dash,
                    "color": line.color,
                },
                label=dict(text=line.label, font=dict(color=line.color)),
            )
            for line in (y_lines or [])
        ]
        + [
            dict(
                type="line",
                yref="paper",
                xref="x",
                x0=line.value,
                y0=0,
                x1=line.value,
                y1=1,
                line={
                    "width": line.width,
                    "dash": line.dash,
                    "color": line.color,
                },
                label=dict(text=line.label, font=dict(color=line.color)),
            )
            for line in (x_lines or [])
        ]
    )

    return LinePlot(**model.__dict__)


def convert_dash_style(cls: type, dash_style: Optional[str]) -> Optional[str]:
    """Convert dash style from Highcharts to Plotly"""
    if dash_style is None:
        return None
    mapping = {
        "Solid": "solid",
        "ShortDash": "dash",
        "ShortDot": "dot",
        "ShortDashDot": "dashdot",
        "ShortDashDotDot": "dashdot",
        "Dot": "dot",
        "Dash": "dash",
        "DashDot": "dashdot",
        "LongDash": "longdash",
        "LongDashDot": "longdashdot",
        "LongDashDotDot": "longdashdot",
    }
    if dash_style in mapping.values():  # Plotly style?
        return dash_style
    elif dash_style in mapping.keys():  # Highcharts style?
        add_validation_warning(
            [LinePlotConfig, cls], f"'{dash_style}' is a deprecated dash style, use '{mapping[dash_style]}'"
        )
        return mapping[dash_style]
    return "solid"


def remove_nones_and_empty_dicts(d: Mapping) -> Dict:
    """Remove None and empty dicts from a dict recursively."""
    return {k: remove_nones_and_empty_dicts(v) for k, v in d.items() if v is not None and v != {}}
