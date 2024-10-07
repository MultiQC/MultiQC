import base64
from functools import lru_cache
import io
import logging
import math
import platform
import queue
import random
import re
from pathlib import Path
import subprocess
from typing import Any, Dict, Generic, List, Mapping, Optional, Tuple, Type, TypeVar, Union

import plotly.graph_objects as go  # type: ignore
from pydantic import BaseModel, ConfigDict, Field, field_serializer, field_validator

from multiqc import config, report
from multiqc.core import tmp_dir
from multiqc.core.strict_helpers import lint_error
from multiqc.plots.plotly import check_plotly_version
from multiqc.types import Anchor, PlotType
from multiqc.utils import mqc_colour
from multiqc.validation import ValidatedConfig, add_validation_warning

logger = logging.getLogger(__name__)

check_plotly_version()


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
    label: Optional[Union[str, Dict[str, Any]]] = None

    @classmethod
    def parse_label(cls, value: Any, _clss: List[Type[Any]]) -> Any:
        if isinstance(value, dict):
            add_validation_warning(
                _clss,
                "Line plot's x_lines or y_lines 'label' field is expected to be a string. "
                "Other fields other than 'text' are deprecated and will be ignored",
            )
            return value["text"]
        return value

    def __init__(self, _clss=None, **data):
        _clss = _clss or []
        if "dashStyle" in data:
            data["dash"] = convert_dash_style(data.pop("dashStyle"), _clss=_clss + [self.__class__])
        if "dash" in data:
            data["dash"] = convert_dash_style(data["dash"], _clss=_clss + [self.__class__])
        super().__init__(**data, _clss=_clss)


class LineBand(ValidatedConfig):
    """
    Extra X1-X2 or Y1-Y2 band added to the plot
    """

    from_: Union[float, int]
    to: Union[float, int]
    color: Optional[str] = None

    def __init__(self, _clss: Optional[List[Type[ValidatedConfig]]] = None, **data: Any):
        super().__init__(**data, _clss=_clss or [self.__class__])


class PConfig(ValidatedConfig):
    id: str
    title: str
    anchor: Optional[Anchor] = None  # unlike id, has to be globally unique
    table_title: Optional[str] = Field(None, deprecated="title")
    height: Optional[int] = None
    width: Optional[int] = None
    square: bool = False
    logswitch: Optional[bool] = False
    logswitch_active: bool = False  # log scale is active
    logswitch_label: str = "Log10"
    cpswitch: Optional[bool] = True
    cpswitch_c_active: bool = True  # percentage scale is _not_ active
    cpswitch_counts_label: str = "Counts"
    cpswitch_percent_label: str = "Percentages"
    xLog: Optional[bool] = Field(None, deprecated="xlog")
    yLog: Optional[bool] = Field(None, deprecated="ylog")
    xlog: bool = False
    ylog: bool = False
    data_labels: List[Union[str, Dict[str, Any]]] = []
    xTitle: Optional[str] = Field(None, deprecated="xlab")
    yTitle: Optional[str] = Field(None, deprecated="ylab")
    xlab: Optional[str] = None
    ylab: Optional[str] = None
    xsuffix: Optional[str] = None
    ysuffix: Optional[str] = None
    tt_suffix: Optional[str] = None
    xLabFormat: Optional[bool] = Field(None, deprecated="xlab_format")
    yLabFormat: Optional[bool] = Field(None, deprecated="ylab_format")
    xLabelFormat: Optional[bool] = Field(None, deprecated="xlab_format")
    yLabelFormat: Optional[bool] = Field(None, deprecated="ylab_format")
    xlab_format: Optional[str] = None
    ylab_format: Optional[str] = None
    tt_label: Optional[str] = None
    xDecimals: Optional[int] = Field(None, deprecated="x_decimals")
    yDecimals: Optional[int] = Field(None, deprecated="y_decimals")
    decimalPlaces: Optional[int] = Field(None, deprecated="tt_decimals")
    x_decimals: Optional[int] = None
    y_decimals: Optional[int] = None
    tt_decimals: Optional[int] = None
    xmin: Optional[Union[float, int]] = None
    xmax: Optional[Union[float, int]] = None
    ymin: Optional[Union[float, int]] = None
    ymax: Optional[Union[float, int]] = None
    xFloor: Optional[Union[float, int]] = Field(None, deprecated="x_clipmin")
    xCeiling: Optional[Union[float, int]] = Field(None, deprecated="x_clipmax")
    yFloor: Optional[Union[float, int]] = Field(None, deprecated="y_clipmin")
    yCeiling: Optional[Union[float, int]] = Field(None, deprecated="y_clipmax")
    x_clipmin: Optional[Union[float, int]] = None
    x_clipmax: Optional[Union[float, int]] = None
    y_clipmin: Optional[Union[float, int]] = None
    y_clipmax: Optional[Union[float, int]] = None
    save_data_file: bool = True
    showlegend: Optional[bool] = None
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
    _actual_cls = None

    @classmethod
    def from_pconfig_dict(cls, pconfig: Union[Mapping[str, Any], "PConfig", None]):
        if pconfig is None:
            lint_error(f"pconfig with required fields 'id' and 'title' must be provided for plot {cls.__name__}")
            return cls(
                id=f"{cls.__name__.lower().replace('config', '')}-{random.randint(1000000, 9999999)}",
                title=cls.__name__,
            )
        elif isinstance(pconfig, PConfig):
            return pconfig
        else:
            return cls(**pconfig)

    def __init__(self, **data: Any):
        super().__init__(**data, _clss=[self.__class__])

        self._actual_cls = self.__class__

        if not self.id:
            self.id = f"{self.__class__.__name__.lower().replace('config', '')}-{random.randint(1000000, 9999999)}"

        # Allow user to overwrite any given config for this plot
        if self.id in config.custom_plot_config:
            for k, v in config.custom_plot_config[self.id].items():
                setattr(self, k, v)

    @classmethod
    def parse_x_bands(cls, data, _clss: List[Type]):
        return [LineBand(**d, _clss=_clss) for d in ([data] if isinstance(data, dict) else data)]

    @classmethod
    def parse_y_bands(cls, data, _clss: List[Type]):
        return [LineBand(**d, _clss=_clss) for d in ([data] if isinstance(data, dict) else data)]

    @classmethod
    def parse_x_lines(cls, data, _clss: List[Type]):
        return [FlatLine(**d, _clss=_clss) for d in ([data] if isinstance(data, dict) else data)]

    @classmethod
    def parse_y_lines(cls, data, _clss: List[Type]):
        return [FlatLine(**d, _clss=_clss) for d in ([data] if isinstance(data, dict) else data)]


class BaseDataset(BaseModel):
    """
    Plot dataset: data and metadata for a single plot. Does not necessarily contain all underlying data,
    as something might be down-sampled for the sake of efficiency of interactive plots. Use intermediate
    data files if you need full data.
    """

    plot_id: str
    label: str
    uid: str
    dconfig: Dict[str, Any]  # user dataset-specific configuration
    layout: Dict[str, Any]  # update when a datasets toggle is clicked, or percentage switch is unselected
    trace_params: Dict[str, Any]
    pct_range: Dict[str, Any]
    n_samples: int

    def create_figure(
        self,
        layout: go.Layout,
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Abstract method to be overridden by specific plots: create a Plotly figure for a dataset, update layout if needed.
        """
        raise NotImplementedError

    def save_data_file(self) -> None:
        """
        Save dataset to disk.
        """
        raise NotImplementedError

    def get_x_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        return None, None

    def get_y_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        return None, None


DatasetT = TypeVar("DatasetT", bound="BaseDataset")
PConfigT = TypeVar("PConfigT", bound="PConfig")


class Plot(BaseModel, Generic[DatasetT, PConfigT]):
    """
    Plot model for serialisation to JSON. Contains enough data to recreate the plot (e.g. in Plotly-JS)
    """

    id: str
    anchor: Anchor  # unlike id, has to be unique
    plot_type: PlotType
    layout: go.Layout
    datasets: List[DatasetT]
    pconfig: PConfigT
    add_log_tab: bool
    add_pct_tab: bool
    l_active: bool
    p_active: bool
    pct_axis_update: Dict[str, Any]
    axis_controlled_by_switches: List[str] = []
    square: bool = False
    flat: bool = False
    defer_render: bool = False

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        use_enum_values=True,
    )

    @field_serializer("layout")
    def serialize_dt(self, layout: go.Layout, _info):
        return layout.to_plotly_json()

    # noinspection PyNestedDecorators
    @field_validator("layout", mode="before")
    @classmethod
    def parse_layout(cls, d: Any):
        if isinstance(d, dict):
            return go.Layout(**d)
        return d

    def __init__(self, **data):
        super().__init__(**data)
        self._set_x_bands_and_range(self.pconfig)
        self._set_y_bands_and_range(self.pconfig)

    @staticmethod
    def initialize(
        plot_type: PlotType,
        pconfig: PConfigT,
        n_samples_per_dataset: List[int],
        id: Optional[str] = None,
        anchor: Optional[str] = None,
        axis_controlled_by_switches: Optional[List[str]] = None,
        default_tt_label: Optional[str] = None,
        defer_render_if_large: bool = True,
        flat_if_very_large: bool = True,
    ) -> "Plot[DatasetT, PConfigT]":
        """
        Initialize a plot model with the given configuration, but without data.
        :param plot_type: plot type
        :param pconfig: plot configuration model
        :param n_samples_per_dataset: number of samples for each dataset, to pre-initialize the base dataset models
        :param id: plot ID
        :param anchor: plot HTML anchor. Unlike ID, must be globally unique
        :param axis_controlled_by_switches: list of axis names that are controlled by the
            log10 scale and percentage switch buttons, e.g. ["yaxis"]
        :param default_tt_label: default tooltip label
        :param defer_render_if_large: whether to defer rendering if the number of data points is large
        :param flat_if_very_large: whether to render flat if the number of data points is very large
        """
        if len(n_samples_per_dataset) == 0:
            raise ValueError("No datasets to plot")

        id = id or pconfig.id
        _anchor: Anchor = Anchor(anchor or pconfig.anchor or id)
        _anchor = Anchor(report.save_htmlid(_anchor))  # make sure it's unique

        # Counts / Percentages / Log10 switch
        add_log_tab: bool = pconfig.logswitch is True and plot_type in [PlotType.BAR, PlotType.LINE]
        add_pct_tab: bool = pconfig.cpswitch is not False and plot_type == PlotType.BAR
        l_active = add_log_tab and pconfig.logswitch_active
        p_active = add_pct_tab and not pconfig.cpswitch_c_active

        height = pconfig.height or 500
        width = pconfig.width
        if pconfig.square:
            width = height

        # Render static image if the number of samples is above the threshold
        flat = False
        if config.plots_force_flat:
            flat = True
        if (
            flat_if_very_large
            and not config.plots_force_interactive
            and n_samples_per_dataset[0] > config.plots_flat_numseries
        ):
            logger.debug(
                f"Plot {id} has {n_samples_per_dataset[0]} samples > config.plots_flat_numseries={config.plots_flat_numseries}, rendering flat"
            )
            flat = True

        defer_render = False
        if defer_render_if_large:
            if (
                n_samples_per_dataset[0] > config.plots_defer_loading_numseries
                or n_samples_per_dataset[0] > config.num_datasets_plot_limit  # DEPRECATED in v1.24
            ):
                logger.debug(
                    f"Plot {id} has {n_samples_per_dataset[0]} samples > config.plots_defer_loading_numseries={config.plots_defer_loading_numseries}, will defer render"
                )
                defer_render = True

        showlegend = pconfig.showlegend
        if showlegend is None:
            showlegend = True if flat else False

        layout: go.Layout = go.Layout(
            title=go.layout.Title(
                text=pconfig.title,
                xanchor="center",
                x=0.5,
                font=dict(size=20),
            ),
            xaxis=go.layout.XAxis(
                gridcolor="rgba(0,0,0,0.05)",
                zerolinecolor="rgba(0,0,0,0.05)",
                color="rgba(0,0,0,0.3)",  # axis labels
                tickfont=dict(size=10, color="rgba(0,0,0,1)"),
                automargin=True,  # auto-expand axis to fit the tick labels
            ),
            yaxis=go.layout.YAxis(
                gridcolor="rgba(0,0,0,0.05)",
                zerolinecolor="rgba(0,0,0,0.05)",
                color="rgba(0,0,0,0.3)",  # axis labels
                tickfont=dict(size=10, color="rgba(0,0,0,1)"),
                automargin=True,  # auto-expand axis to fit the tick labels
            ),
            height=height,
            width=width,
            paper_bgcolor="white",
            plot_bgcolor="white",
            font=dict(family="'Lucida Grande', 'Open Sans', verdana, arial, sans-serif"),
            colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
            autosize=True,
            margin=go.layout.Margin(
                pad=5,  # pad sample names in a bar graph a bit
                t=50,  # more compact title
                r=15,  # remove excessive whitespace on the right
                b=65,  # remove excessive whitespace on the bottom
                l=60,  # remove excessive whitespace on the left
            ),
            hoverlabel=go.layout.Hoverlabel(
                namelength=-1,  # do not crop sample names inside hover label <extra></extra>
            ),
            modebar=go.layout.Modebar(
                bgcolor="rgba(0, 0, 0, 0)",
                color="rgba(0, 0, 0, 0.5)",
                activecolor="rgba(0, 0, 0, 1)",
            ),
            showlegend=showlegend,
            legend=go.layout.Legend(
                orientation="h",
                yanchor="top",
                y=-0.15,
                xanchor="center",
                x=0.5,
            )
            if flat
            else None,
        )
        # Layout update for the counts/percentage switch
        pct_axis_update = dict(
            ticksuffix="%",
            hoverformat=".1f",
        )

        axis_controlled_by_switches = axis_controlled_by_switches or []
        if pconfig.xlog:
            _set_axis_log_scale(layout.xaxis)
            if "xaxis" in axis_controlled_by_switches:
                axis_controlled_by_switches.remove("xaxis")
        if pconfig.ylog:
            _set_axis_log_scale(layout.yaxis)
            if "yaxis" in axis_controlled_by_switches:
                axis_controlled_by_switches.remove("yaxis")
        if add_log_tab and l_active:
            for axis in axis_controlled_by_switches:
                layout[axis].type = "log"

        datasets = []
        for idx, n_samples in enumerate(n_samples_per_dataset):
            dataset = BaseDataset(
                plot_id=id,
                label=str(idx + 1),
                uid=id,
                dconfig=dict(),
                layout=dict(),
                trace_params=dict(),
                pct_range=dict(  # range for the percentage view for each axis
                    xaxis=dict(min=0, max=100),
                    yaxis=dict(min=0, max=100),
                ),
                n_samples=n_samples,
            )
            if len(n_samples_per_dataset) > 1:
                dataset.uid += f"_{idx + 1}"

            data_label: Union[str, Dict[str, Union[str, Dict[str, str]]]] = (
                pconfig.data_labels[idx] if idx < len(pconfig.data_labels) else {}
            )
            dconfig: Dict[str, Union[str, Dict[str, str]]] = (
                data_label if isinstance(data_label, dict) else {"name": data_label}
            )
            label = dconfig.get("name", dconfig.get("label", str(idx + 1)))
            assert isinstance(label, str)
            dataset.label = label
            if "ylab" not in dconfig and not pconfig.ylab:
                dconfig["ylab"] = dconfig.get("name", dconfig.get("label", ""))

            if "title" not in dconfig:
                dconfig["title"] = pconfig.title
            subtitles = []
            if len(n_samples_per_dataset) > 1:
                subtitles += [dataset.label]
            if n_samples > 1:
                subtitles += [f"{n_samples} samples"]
            if subtitles:
                title = dconfig.get("title", "")
                assert isinstance(title, str)
                title += f"<br><sup>{', '.join(subtitles)}</sup>"
                dconfig["title"] = title

            dataset.layout, dataset.trace_params = _dataset_layout(pconfig, dconfig, default_tt_label)
            dataset.dconfig = dconfig
            datasets.append(dataset)

        return Plot(
            plot_type=plot_type,
            pconfig=pconfig,
            id=id,
            anchor=_anchor,
            datasets=datasets,
            layout=layout,
            add_log_tab=add_log_tab,
            add_pct_tab=add_pct_tab,
            l_active=l_active,
            p_active=p_active,
            pct_axis_update=pct_axis_update,
            axis_controlled_by_switches=axis_controlled_by_switches,
            square=pconfig.square,
            flat=flat,
            defer_render=defer_render,
        )

    def _set_x_bands_and_range(self, pconfig: PConfigT):
        x_minrange = pconfig.x_minrange
        x_bands = pconfig.x_bands
        x_lines = pconfig.x_lines

        if x_bands or x_lines or x_minrange:
            # same as above but for x-axis
            for dataset in self.datasets:
                minval = dataset.layout["xaxis"]["autorangeoptions"]["minallowed"]
                maxval = dataset.layout["xaxis"]["autorangeoptions"]["maxallowed"]
                dminval, dmaxval = dataset.get_x_range()
                if dminval is not None:
                    minval = min(minval, dminval) if minval is not None else dminval
                if dmaxval is not None:
                    maxval = max(maxval, dmaxval) if maxval is not None else dmaxval
                clipmin = dataset.layout["xaxis"]["autorangeoptions"]["clipmin"]
                clipmax = dataset.layout["xaxis"]["autorangeoptions"]["clipmax"]
                if clipmin is not None and minval is not None and clipmin > minval:
                    minval = clipmin
                if clipmax is not None and maxval is not None and clipmax < maxval:
                    maxval = clipmax
                if x_minrange is not None and maxval is not None and minval is not None:
                    maxval = max(maxval, minval + x_minrange)
                if self.layout.xaxis.type == "log":
                    minval = math.log10(minval) if minval is not None and minval > 0 else None
                    maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
                dataset.layout["xaxis"]["range"] = [minval, maxval]

        if not self.layout.shapes:
            self.layout.shapes = []
        self.layout.shapes = (
            list(self.layout.shapes)
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

    def _set_y_bands_and_range(self, pconfig: PConfigT):
        y_minrange = pconfig.y_minrange
        y_bands = pconfig.y_bands
        y_lines = pconfig.y_lines

        if y_bands or y_lines or y_minrange:
            # We don't want the bands to affect the calculated axis range, so we
            # find the min and the max from data points, and manually set the range.
            for dataset in self.datasets:
                minval = dataset.layout["yaxis"]["autorangeoptions"]["minallowed"]
                maxval = dataset.layout["yaxis"]["autorangeoptions"]["maxallowed"]
                dminval, dmaxval = dataset.get_y_range()
                if dminval is not None:
                    minval = min(minval, dminval) if minval is not None else dminval
                if dmaxval is not None:
                    maxval = max(maxval, dmaxval) if maxval is not None else dmaxval
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
                if self.layout.yaxis.type == "log":
                    minval = math.log10(minval) if minval is not None and minval > 0 else None
                    maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
                dataset.layout["yaxis"]["range"] = [minval, maxval]

        if not self.layout.shapes:
            self.layout.shapes = []
        self.layout.shapes = (
            list(self.layout.shapes)
            + [
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
        )

    def show(self, dataset_id: Union[int, str] = 0, flat: bool = False, **kwargs):
        """
        Show the plot in an interactive environment such as Jupyter notebook.

        @param dataset_id: index of the dataset to plot
        @param flat: whether to save a flat image or an interactive plot
        """
        fig = self.get_figure(dataset_id=dataset_id, flat=flat, **kwargs)
        if flat:
            try:
                from IPython.core.display import HTML  # type: ignore
            except ImportError:
                raise ImportError(
                    "IPython is required to show plot. The function is expected to be run in an interactive environment, "
                    "such as Jupyter notebook. To save plot to file, use Plot.save method"
                )

            return HTML(
                fig_to_static_html(
                    fig,
                    active=True,
                    embed_in_html=True,
                    export_plots=False,
                    file_name=self.id,
                )
            )
        else:
            return fig

    @staticmethod
    def _proc_save_args(filename: str, flat: Optional[bool]) -> Tuple[str, bool]:
        if isinstance(filename, (Path, str)):
            if Path(filename).suffix.lower() == ".html":
                if flat is not None and flat is True:
                    raise ValueError("Set flat=False to save an interactive plot as an HTML file")
                flat = False
            else:
                if flat is not None and flat is False:
                    raise ValueError("Set flat=True to save a static plot as an image file")
                flat = True
        return filename, flat

    def save(self, filename, dataset_id: Union[int, str] = 0, flat=None, **kwargs):
        """
        Save the plot to a file. Will write an HTML with an interactive plot -
        unless flat=True is specified, in which case will write a PNG file.

        @param filename: a string representing a local file path or a writeable object
        (e.g. a pathlib.Path object or an open file descriptor). If the filename ends with ".html",
        an interactive plot will be saved, otherwise a flat image.
        @param dataset_id: index of the dataset to plot
        @param flat: whether to save a static image instead of an interactive HTML.
        """
        filename, flat = self._proc_save_args(filename, flat)

        fig = self.get_figure(dataset_id=dataset_id, flat=flat, **kwargs)
        if flat:
            fig.write_image(
                filename,
                scale=2,
                width=fig.layout.width,
                height=fig.layout.height,
            )
        else:
            fig.write_html(
                filename,
                include_plotlyjs="cdn",
                full_html=False,
            )
        logger.info(f"Plot saved to {filename}")

    def get_figure(
        self,
        dataset_id: Union[int, str],
        is_log: bool = False,
        is_pct: bool = False,
        flat: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Public method: create a Plotly Figure object.
        """
        if isinstance(dataset_id, str):
            for i, d in enumerate(self.datasets):
                if d.label == dataset_id:
                    dataset = d
                    break
            else:
                dataset = self.datasets[0]
        else:
            dataset = self.datasets[dataset_id]

        layout = go.Layout(self.layout.to_plotly_json())  # make a copy
        layout.update(**dataset.layout)
        if flat:
            if config.simple_output:
                layout.width = 600
            else:
                layout.width = 1100
        for axis in self.axis_controlled_by_switches:
            layout[axis].type = "linear"
            minval = layout[axis].autorangeoptions["minallowed"]
            maxval = layout[axis].autorangeoptions["maxallowed"]
            if is_pct:
                layout[axis].update(self.pct_axis_update)
                minval = dataset.pct_range.get(axis, {}).get("min", 0)
                maxval = dataset.pct_range.get(axis, {}).get("max", 100)
            if is_log:
                layout[axis].type = "log"
                minval = math.log10(minval) if minval is not None and minval > 0 else None
                maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
            layout[axis].autorangeoptions["minallowed"] = minval
            layout[axis].autorangeoptions["maxallowed"] = maxval
        return dataset.create_figure(layout, is_log, is_pct, **kwargs)

    def __repr__(self):
        d = {k: v for k, v in self.__dict__.items() if k not in ("datasets", "layout")}
        return f"<{self.__class__.__name__} {self.id} {d}>"

    def add_to_report(self, plots_dir_name: Optional[str] = None) -> str:
        """
        Build and add the plot data to the report, return an HTML wrapper.
        """
        for ds in self.datasets:
            ds.uid = self.id
            if len(self.datasets) > 1:  # for flat plots, each dataset will have its own unique ID
                ds.uid = report.save_htmlid(f"{self.id}_{ds.label}", skiplint=True)

        if self.flat:
            if is_running_under_rosetta():
                # Kaleido is unstable under rosetata, falling back to interactive plots
                return self.interactive_plot()
            try:
                html = self.flat_plot(plots_dir_name=plots_dir_name)
            except ValueError:
                logger.error(f"Unable to export plot to flat images: {self.id}, falling back to interactive plot")
                html = self.interactive_plot()
        else:
            html = self.interactive_plot()
            if config.export_plots and not is_running_under_rosetta():
                try:
                    self.flat_plot(embed_in_html=False, plots_dir_name=plots_dir_name)
                except ValueError:
                    logger.error(f"Unable to export plot to flat images: {self.id}")

        return html

    def save_data_files(self):
        if self.id != "general_stats_table" and self.pconfig.save_data_file:
            for dataset in self.datasets:
                dataset.save_data_file()

    def interactive_plot(self) -> str:
        html = '<div class="mqc_hcplot_plotgroup">'

        html += self.__control_panel(flat=False)

        # This width only affects the space before plot is rendered, and the initial
        # height for the resizing function. For the actual plot container, Plotly will
        # re-calculate the wrapper size after rendering.
        height_style = f'style="height:{self.layout.height + 7}px"' if self.layout.height else ""
        defer_render_style = "defer_render" if self.defer_render else ""
        html += f"""
        <div class="hc-plot-wrapper hc-{self.plot_type}-wrapper" id="{self.anchor}-wrapper" {height_style}>
            <div
                id="{self.anchor}"
                class="hc-plot hc-{self.plot_type}-plot not_loaded not_rendered {defer_render_style}">
            </div>
            <div class="created-with-multiqc">Created with MultiQC</div>
        </div>"""
        html += "</div>"

        # Saving compressed data for JavaScript to pick up and uncompress.
        report.plot_data[self.anchor] = self.model_dump(warnings=False)
        return html

    def flat_plot(self, embed_in_html: Optional[bool] = None, plots_dir_name: Optional[str] = None) -> str:
        embed_in_html = embed_in_html if embed_in_html is not None else not config.development
        if not embed_in_html and plots_dir_name is None:
            raise ValueError("plots_dir_name is required for non-embedded plots")

        html = "".join(
            [
                '<p class="text-info">',
                "<small>" '<span class="glyphicon glyphicon-picture" aria-hidden="true"></span> ',
                "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work ",
                '(see the <a href="https://docs.seqera.io/multiqc/development/plots/#interactive--flat-image-plots" target="_blank">docs</a>).',
                "</small>",
                "</p>",
            ]
        )
        html += f'<div class="mqc_mplplot_plotgroup" id="plotgroup-{self.anchor}" data-plot-anchor={self.anchor}>'

        if not config.simple_output:
            html += self.__control_panel(flat=True)

        # Go through datasets creating plots
        for ds_idx, dataset in enumerate(self.datasets):
            html += fig_to_static_html(
                self.get_figure(ds_idx, flat=True),
                active=ds_idx == 0 and not self.p_active and not self.l_active,
                file_name=dataset.uid if not self.add_log_tab and not self.add_pct_tab else f"{dataset.uid}-cnt",
                plots_dir_name=plots_dir_name,
                embed_in_html=embed_in_html,
            )
            if self.add_pct_tab:
                html += fig_to_static_html(
                    self.get_figure(ds_idx, is_pct=True, flat=True),
                    active=ds_idx == 0 and self.p_active,
                    file_name=f"{dataset.uid}-pct",
                    plots_dir_name=plots_dir_name,
                    embed_in_html=embed_in_html,
                )
            if self.add_log_tab:
                html += fig_to_static_html(
                    self.get_figure(ds_idx, is_log=True, flat=True),
                    active=ds_idx == 0 and self.l_active,
                    file_name=f"{dataset.uid}-log",
                    plots_dir_name=plots_dir_name,
                    embed_in_html=embed_in_html,
                )
            if self.add_pct_tab and self.add_log_tab:
                html += fig_to_static_html(
                    self.get_figure(ds_idx, is_pct=True, is_log=True, flat=True),
                    active=ds_idx == 0 and self.p_active and self.l_active,
                    file_name=f"{dataset.uid}-pct-log",
                    plots_dir_name=plots_dir_name,
                    embed_in_html=embed_in_html,
                )

        html += "</div>"
        return html

    def _btn(
        self,
        cls: str,
        label: str,
        data_attrs: Optional[Dict[str, str]] = None,
        attrs: Optional[Dict[str, str]] = None,
        pressed: bool = False,
    ) -> str:
        """
        Build a switch button for the plot.
        """
        attrs = attrs.copy() if attrs else {}
        attrs_str = " ".join([f'{k}="{v}"' for k, v in attrs.items()])

        data_attrs = data_attrs.copy() if data_attrs else {}
        if "plot-anchor" not in data_attrs:
            data_attrs["plot-anchor"] = self.anchor
        data_attrs_str = " ".join([f'data-{k}="{v}"' for k, v in data_attrs.items()])

        return f'<button {attrs_str} class="btn btn-default btn-sm {cls} {"active" if pressed else ""}" {data_attrs_str}>{label}</button>\n'

    def buttons(self, flat: bool) -> List[str]:
        """
        Build buttons for control panel
        """
        switch_buttons = ""
        cls = "mpl_switch_group" if flat else "interactive-switch-group"
        # Counts / percentages / log10 switches
        if self.add_pct_tab or self.add_log_tab:
            if self.add_pct_tab:
                switch_buttons += self._btn(
                    cls=f"{cls} percent-switch",
                    label=self.pconfig.cpswitch_percent_label,
                    pressed=self.p_active,
                )
            if self.add_log_tab:
                switch_buttons += self._btn(
                    cls=f"{cls} log10-switch",
                    label=self.pconfig.logswitch_label,
                    pressed=self.l_active,
                )

        # Buttons to cycle through different datasets
        if len(self.datasets) > 1:
            switch_buttons += f'<div class="btn-group {cls} dataset-switch-group">\n'
            for ds_idx, ds in enumerate(self.datasets):
                data_attrs: Dict[str, str] = {
                    "dataset-index": str(ds_idx),
                    # For flat plots, we will generate separate flat images for each
                    # dataset and view, so have to save individual image IDs.
                    "dataset-uid": ds.uid,
                }
                switch_buttons += self._btn(
                    cls="mr-auto",
                    label=ds.label,
                    data_attrs=data_attrs,
                    pressed=ds_idx == 0,
                )
            switch_buttons += "</div>\n\n"

        export_btn = ""
        if not flat:
            export_btn = self._btn(
                cls="export-plot",
                label="Export Plot",
                data_attrs={"plot-anchor": str(self.anchor), "type": str(self.plot_type)},
            )
        return [switch_buttons, export_btn]

    def __control_panel(self, flat: bool) -> str:
        """
        Add buttons: percentage on/off, log scale on/off, datasets switch panel
        """
        buttons = "\n".join(self.buttons(flat=flat))
        html = f"<div class='row'>\n<div class='col-xs-12'>\n{buttons}\n</div>\n</div>\n\n"
        return html


def _export_plot(fig, file_ext, plot_path, write_kwargs) -> Optional[str]:
    if is_running_under_rosetta():
        return None

    try:
        if file_ext == "svg":
            # Cannot add logo to SVGs
            fig.write_image(plot_path, **write_kwargs)
        else:
            img_buffer = io.BytesIO()
            fig.write_image(img_buffer, **write_kwargs)
            img_buffer = add_logo(img_buffer, format=file_ext)
            with open(plot_path, "wb") as f:
                f.write(img_buffer.getvalue())
            img_buffer.close()
    except Exception as e:
        logger.error(f"Unable to export plot to {file_ext.upper()} image: {e}")
        return None
    else:
        return plot_path


def _export_plot_to_buffer(fig, write_kwargs) -> Optional[str]:
    try:
        img_buffer = io.BytesIO()
        fig.write_image(img_buffer, **write_kwargs)
        img_buffer = add_logo(img_buffer, format="PNG")
        # Convert to a base64 encoded string
        b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
        img_src = f"data:image/png;base64,{b64_img}"
        img_buffer.close()
    except Exception as e:
        logger.error(f"Unable to export PNG figure to static image: {e}")
        return None
    else:
        return img_src


plot_export_has_failed: bool = False


def fig_to_static_html(
    fig: go.Figure,
    active: bool = True,
    export_plots: Optional[bool] = None,
    embed_in_html: Optional[bool] = None,
    plots_dir_name: Optional[str] = None,
    file_name: Optional[str] = None,
) -> str:
    """
    Build one static image, return an HTML wrapper.
    """
    if is_running_under_rosetta():
        raise ValueError(
            "Detected Rosetta process, meaning running in an x86_64 container hosted by Apple Silicon. "
            "Plot export is unstable and will be skipped"
        )

    global plot_export_has_failed
    if plot_export_has_failed:
        raise ValueError("Could not previously export a plots, so won't try again")

    embed_in_html = embed_in_html if embed_in_html is not None else not config.development
    export_plots = export_plots if export_plots is not None else config.export_plots

    assert fig.layout.width
    scale = 2.0  # higher detail (to look sharp on the retina display)
    scale *= config.plots_export_font_scale  # bigger font if configured in the settings
    write_kwargs = dict(
        width=fig.layout.width / config.plots_export_font_scale,  # While interactive plots take full width of screen,
        # for the flat plots we explicitly set width
        height=fig.layout.height / config.plots_export_font_scale,
        scale=scale,  # higher detail (retina display)
        # engine="orca",  # kaleido gets frozen in docker environments
    )

    formats = set(config.export_plot_formats) if export_plots else set()
    if not embed_in_html and "png" not in formats:
        if not export_plots:
            formats = {"png"}
        else:
            formats.add("png")

    # Save the plot to the data directory if export is requested
    png_is_written = False

    if formats:
        if file_name is None:
            raise ValueError("file_name is required for export_plots")
        for file_ext in formats:
            plot_path = tmp_dir.plots_tmp_dir() / file_ext / f"{file_name}.{file_ext}"
            plot_path.parent.mkdir(parents=True, exist_ok=True)

            try:
                if not plot_export_has_failed:
                    # Running for the first time, so doing a safe run in a subprocess to find out if it freezes the process or now
                    _export_plot(fig, file_ext, plot_path, write_kwargs)
            except Exception as e:
                msg = f"{file_name}: Unable to export plot to {file_ext.upper()} image"
                logger.error(f"{msg}. {e}")
                plot_export_has_failed = True
                raise ValueError(msg)  # Raising to the caller to fall back to interactive plots
            else:
                if file_ext == "png":
                    png_is_written = True

    # Now writing the PNGs for the HTML
    if not embed_in_html:
        if file_name is None:
            raise ValueError("file_name is required for non-embedded plots")
        if plots_dir_name is None:
            raise ValueError("plots_dir_name is required for non-embedded plots")
        # Using file written in the config.export_plots block above
        img_path = Path(plots_dir_name) / "png" / f"{file_name}.png"
        if not png_is_written:  # Could not write in the block above
            raise ValueError(f"Unable to export plot to PNG image {file_name}")
        img_src = str(img_path)
    else:
        _img_src = _export_plot_to_buffer(fig, write_kwargs)
        if _img_src is None:
            raise ValueError("Unable to export PNG figure to static image")
        img_src = _img_src

    # Should this plot be hidden on report load?
    style = "" if active else "display:none;"
    return "".join(
        [
            f'<div class="mqc_mplplot" style="{style}" id="{file_name}">',
            f'<img src="{img_src}" height="{fig.layout.height}px" width="{fig.layout.width}px"/>',
            "</div>",
        ]
    )


def add_logo(
    img_buffer: io.BytesIO,
    format: str = "png",
    text: str = "Created with MultiQC",
    font_size: int = 16,
) -> io.BytesIO:
    try:
        from PIL import Image, ImageDraw

        # Load the image from the BytesIO object
        image = Image.open(img_buffer)

        # Create a drawing context
        draw = ImageDraw.Draw(image)

        # Define the text position. In order to do that, first calculate the expected
        # text block width, given the font size.
        # noinspection PyArgumentList
        text_width: float = draw.textlength(text, font_size=font_size)
        position: Tuple[int, int] = (image.width - int(text_width) - 3, image.height - 30)

        # Draw the text
        draw.text(position, text, fill="#9f9f9f", font_size=font_size)

        # Save the image to a BytesIO object
        output_buffer = io.BytesIO()
        image.save(output_buffer, format=format)
        output_buffer.seek(0)

    except Exception as e:
        logger.warning(f"Failure adding logo to the plot: {e}")
        output_buffer = img_buffer

    return output_buffer


def _set_axis_log_scale(axis):
    axis.type = "log"
    minval = axis.autorangeoptions["minallowed"]
    maxval = axis.autorangeoptions["maxallowed"]
    minval = math.log10(minval) if minval is not None and minval > 0 else None
    maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
    axis.autorangeoptions["minallowed"] = minval
    axis.autorangeoptions["maxallowed"] = maxval


def rename_deprecated_highcharts_keys(conf: Dict) -> Dict:
    """
    Rename the deprecated HighCharts-specific terminology in a config.
    """
    conf = conf.copy()
    if "yCeiling" in conf:
        conf["yaxis"] = conf.pop("y_ceiling")
    if "xAxis" in conf:
        conf["xaxis"] = conf.pop("xAxis")
    if "tooltip" in conf:
        conf["hovertemplate"] = conf.pop("tooltip")
    return conf


def _dataset_layout(
    pconfig: PConfig,
    dconfig: Dict,
    default_tt_label: Optional[str] = None,
) -> Tuple[Dict, Dict]:
    """
    Given plot config and dataset config, set layout and trace params.
    """
    pconfig = pconfig.model_copy()
    for k, v in dconfig.items():
        if k in pconfig.model_fields:
            setattr(pconfig, k, v)

    ysuffix = pconfig.ysuffix if pconfig.ysuffix is not None else pconfig.tt_suffix
    xsuffix = pconfig.xsuffix

    # Handle hover tooltip options deprecated in 1.21:
    if ysuffix is None and pconfig.ylab_format:
        if "}" in pconfig.ylab_format:
            ysuffix = pconfig.ylab_format.split("}")[1]
    if xsuffix is None and pconfig.xlab_format:
        if "}" in pconfig.xlab_format:
            xsuffix = pconfig.xlab_format.split("}")[1]

    # Set or remove space in known suffixes
    KNOWN_SUFFIXES = ["%", "x", "X", "k", "M", " bp", " kbp", " Mbp"]
    for suf in KNOWN_SUFFIXES:
        if ysuffix is not None and ysuffix == suf.strip():
            ysuffix = suf
        if xsuffix is not None and xsuffix == suf.strip():
            xsuffix = suf

    # Set % suffix from ylab if it's in form like "% reads"
    if ysuffix is None and pconfig.ylab:
        if "%" in pconfig.ylab or "percentage" in pconfig.ylab.lower():
            ysuffix = "%"
        for suf in KNOWN_SUFFIXES:
            if pconfig.ylab.endswith(f" ({suf.strip()})"):
                ysuffix = suf
    if xsuffix is None and pconfig.xlab:
        if "%" in pconfig.xlab or "percentage" in pconfig.xlab.lower():
            xsuffix = "%"
        for suf in KNOWN_SUFFIXES:
            if pconfig.xlab.endswith(f" ({suf.strip()})"):
                xsuffix = suf

    tt_label: Optional[str] = None
    if pconfig.tt_label is not None:
        # Clean the hover tooltip label, add missing <br> into the beginning, populate suffixes if missing
        tt_label = pconfig.tt_label
        tt_label = _clean_config_tt_label(tt_label)

        if ysuffix is None or xsuffix is None:
            # if "%" or other suffix is in the hover label, parse that suffix to add it to the ticks
            parts = tt_label.split("%{")
            for part in parts:
                if ysuffix is None and part.startswith("y") and "}" in part:
                    info = part.split("}")[1].replace("</b>", "")
                    info = info.split(":")[0].split(",")[0].strip().split(" ")[0]
                    if info:
                        for suf in KNOWN_SUFFIXES:
                            if info == suf.strip():
                                ysuffix = suf
                                break
                elif xsuffix is None and part.startswith("x") and "}" in part:
                    info = part.split("}")[1].replace("</b>", "")
                    info = info.split(":")[0].split(",")[0].strip().split(" ")[0]
                    if info:
                        for suf in KNOWN_SUFFIXES:
                            if info == suf.strip():
                                xsuffix = suf
                                break

        # As the suffix will be added automatically for the simple format ({y}), remove it from the label
        if ysuffix is not None:
            if "{y}" + ysuffix in tt_label:
                tt_label = tt_label.replace("{y}" + ysuffix, "{y}")
            if "{y} " + ysuffix in tt_label:
                tt_label = tt_label.replace("{y} " + ysuffix, "{y}")
        if xsuffix is not None:
            if "{x}" + xsuffix in tt_label:
                tt_label = tt_label.replace("{x}" + xsuffix, "{x}")
            if "{x} " + xsuffix in tt_label:
                tt_label = tt_label.replace("{x} " + xsuffix, "{x}")

        # add missing line break between the sample name and the key-value pair
        if not tt_label.startswith("<br>"):
            tt_label = "<br>" + tt_label
    elif default_tt_label is not None:
        tt_label = default_tt_label

    if tt_label:
        hovertemplate = "<b>%{text}</b>" + tt_label + "<extra></extra>"
    else:
        hovertemplate = None

    # `hoverformat` describes how plain "{y}" or "{x}" are formatted in `hovertemplate`
    y_decimals = pconfig.tt_decimals if pconfig.tt_decimals is not None else pconfig.y_decimals
    y_hoverformat = f",.{y_decimals}f" if y_decimals is not None else None

    x_decimals = pconfig.x_decimals
    x_hoverformat = f",.{x_decimals}f" if x_decimals is not None else None

    layout = dict(
        title=dict(text=pconfig.title),
        xaxis=dict(
            hoverformat=x_hoverformat,
            ticksuffix=xsuffix or "",
            title=dict(text=pconfig.xlab),
            rangemode="tozero" if pconfig.xmin == 0 else "normal",
            autorangeoptions=dict(
                clipmin=pconfig.x_clipmin,
                clipmax=pconfig.x_clipmax,
                minallowed=pconfig.xmin,
                maxallowed=pconfig.xmax,
            ),
        ),
        yaxis=dict(
            hoverformat=y_hoverformat,
            ticksuffix=ysuffix or "",
            title=dict(text=pconfig.ylab),
            rangemode="tozero" if pconfig.ymin == 0 else "normal",
            autorangeoptions=dict(
                clipmin=pconfig.y_clipmin,
                clipmax=pconfig.y_clipmax,
                minallowed=pconfig.ymin,
                maxallowed=pconfig.ymax,
            ),
        ),
    )

    trace_params = {}
    if hovertemplate:
        trace_params["hovertemplate"] = hovertemplate

    return layout, trace_params


def _clean_config_tt_label(tt_label: str) -> str:
    replace_d = {  # Convert HighCharts format to Plotly format
        "{point.x": "%{x",
        "{point.y": "%{y",
        "x:>": "x:",
        "y:>": "y:",
        "{point.category}": "%{x}",
        "<strong>": "<b>",
        "</strong>": "</b>",
        "<br/>": "<br>",
    }
    for k, v in replace_d.items():
        tt_label = tt_label.replace(k, v)
    return tt_label


def split_long_string(s: str, max_width=80) -> List[str]:
    """
    Split string into lines of max_width characters
    """
    lines = []
    current_line = ""
    words = re.split(r"(\W+)", s)
    for word in words:
        if len(current_line + word) <= max_width:
            current_line += word
        else:
            if current_line:
                lines.append(current_line)
            current_line = word

    if current_line:
        lines.append(current_line)

    return lines


def convert_dash_style(dash_style: Optional[str], _clss: Optional[List[type]] = None) -> Optional[str]:
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
        add_validation_warning(_clss or [], f"'{dash_style}' is a deprecated dash style, use '{mapping[dash_style]}'")
        return mapping[dash_style]
    return "solid"


ROSETTA_WARNING = (
    "\n===========================\n"
    "WARNING: Detected a Rosetta process, implying that MultiQC is running inside an x86_64 container, but hosted "
    "by Apple Silicon.\n\nStatic plot export uses Kaleido - a tool that is unstable under conflicting architectures, so "
    "MultiQC will not save static PNG/PDF/SVG plots to disk. If you still need those, please use a container with "
    "a compatible architecture. The official Dockerhub image multiqc/multiqc:latest is available for both "
    "platforms, and can be pulled with:\n\n"
    "docker pull multiqc/multiqc:latest\n"
    "==========================="
)


@lru_cache()
def is_running_under_rosetta() -> bool:
    """
    Detect if running in an x86_64 container hosted by Apple Silicon. Kaleido often freezes in such containers:
    https://github.com/MultiQC/MultiQC/issues/2667
    https://github.com/MultiQC/MultiQC/issues/2812
    https://github.com/MultiQC/MultiQC/issues/2867
    So we have to know when to disable plot export and show a warning.
    """
    # If not x86_64 architecture, rosetta wouldn't be needed
    if platform.machine() != "x86_64":
        return False

    # Check for Rosetta processes
    try:
        output = subprocess.check_output(["ps", "aux"], universal_newlines=True)
        if "/run/rosetta/rosetta" in output:
            logger.warning(ROSETTA_WARNING)
            return True
    except subprocess.CalledProcessError:
        pass

    return False
