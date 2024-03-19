import base64
import dataclasses
import io
import logging
import random
import re
import string
from abc import abstractmethod, ABC
from enum import Enum
from pathlib import Path
from typing import Dict, Union, List, Optional, Tuple

import math
import plotly.graph_objects as go

from multiqc.plots.plotly import check_plotly_version
from multiqc.utils import mqc_colour, config, report

logger = logging.getLogger(__name__)


class PlotType(Enum):
    """
    Plot labels used in custom content, as well as in JS to load plot into Plotly-js
    """

    LINE = "xy_line"
    VIOLIN = "violin"
    BAR = "bar_graph"
    SCATTER = "scatter"
    HEATMAP = "heatmap"
    BOX = "box"


@dataclasses.dataclass
class BaseDataset(ABC):
    """
    Structured version of dataset config dictionary
    """

    plot_id: str
    label: str
    uid: str
    dconfig: Dict  # user dataset-specific configuration
    layout: Dict  # update when a datasets toggle is clicked, or percentage switch is unselected
    trace_params: Dict
    pct_range: Dict

    def dump_for_javascript(self) -> Dict:
        return {k: v for k, v in self.__dict__.items()}

    def create_figure(self, layout: go.Layout, is_log=False, is_pct=False):
        """
        To be overridden by specific plots: create a Plotly figure for a dataset, update layout if needed.
        """
        raise NotImplementedError


class Plot(ABC):
    """Structured version of plot config dictionary"""

    # Default width for flat plots
    FLAT_PLOT_WIDTH = 1100

    def __init__(
        self,
        plot_type: PlotType,
        pconfig: Dict,
        n_datasets: int,
        id: Optional[str] = None,
    ):
        check_plotly_version()

        if n_datasets == 0:
            raise ValueError("No datasets to plot")

        self.id = id or pconfig.get("id")
        if self.id is None:  # ID of the plot group
            uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
            is_static_suf = "static_" if config.plots_force_flat else ""
            self.id = f"mqc_{is_static_suf}plot_{uniq_suffix}"

        self.plot_type = plot_type
        self.flat = config.plots_force_flat or (
            not config.plots_force_interactive and n_datasets > config.plots_flat_numseries
        )
        self.pconfig = pconfig

        # Counts / Percentages / Log10 switch
        self.add_log_tab = pconfig.get("logswitch", False) and plot_type in [PlotType.BAR, PlotType.LINE]
        self.add_pct_tab = pconfig.get("cpswitch", True) is not False and plot_type == PlotType.BAR
        self.l_active = self.add_log_tab and pconfig.get("logswitch_active") is True
        self.p_active = self.add_pct_tab and pconfig.get("cpswitch_c_active", True) is not True

        # Per-dataset configurations
        self.datasets: List[BaseDataset] = [
            BaseDataset(
                plot_id=self.id,
                label=str(i + 1),
                uid=self.id,
                dconfig=dict(),
                layout=dict(),
                trace_params=dict(),
                pct_range=dict(  # range for the percentage view for each axis
                    xaxis=dict(min=0, max=100),
                    yaxis=dict(min=0, max=100),
                ),
            )
            for i in range(n_datasets)
        ]

        height = self.pconfig.get("height", 500)
        width = self.pconfig.get("width")
        if self.pconfig.get("square"):
            width = height

        title = self.pconfig.get("table_title", self.pconfig.get("title"))
        if not title:
            if config.strict:
                errmsg = f"LINT: 'title' is missing from plot config for plot '{self.id}'"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)

        self.layout = go.Layout(
            title=go.layout.Title(
                text=title,
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
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
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
            showlegend=self.flat,
            legend=go.layout.Legend(
                orientation="h",
                yanchor="top",
                y=-0.15,
                xanchor="center",
                x=0.5,
            )
            if self.flat
            else None,
        )
        # Layout update for the counts/percentage switch
        self.pct_axis_update = dict(
            ticksuffix="%",
            hoverformat=".1f",
        )

        self._axis_controlled_by_switches = self.axis_controlled_by_switches()
        if self.pconfig.get("xlog", self.pconfig.get("xLog")):
            self._set_axis_log_scale(self.layout.xaxis)
            if "xaxis" in self._axis_controlled_by_switches:
                self._axis_controlled_by_switches.remove("xaxis")
        if self.pconfig.get("ylog", self.pconfig.get("yLog")):
            self._set_axis_log_scale(self.layout.yaxis)
            if "yaxis" in self._axis_controlled_by_switches:
                self._axis_controlled_by_switches.remove("yaxis")
        if self.add_log_tab and self.l_active:
            for axis in self._axis_controlled_by_switches:
                self.layout[axis].type = "log"

        dconfigs: List[Union[str, Dict[str, str]]] = pconfig.get("data_labels", [])
        for idx, dataset in enumerate(self.datasets):
            if n_datasets > 1:
                dataset.uid += f"_{idx + 1}"

            if idx < len(dconfigs):
                dconfig = dconfigs[idx]
            else:
                dconfig = {}

            if not isinstance(dconfig, (str, dict)):
                logger.warning(f"Invalid data_labels type: {type(dconfig)}. Must be a string or a dict.")
            dconfig = dconfig if isinstance(dconfig, dict) else {"name": dconfig}
            dataset.label = dconfig.get("name", dconfig.get("label", idx + 1))
            if "ylab" not in dconfig and "ylab" not in self.pconfig:
                dconfig["ylab"] = dconfig.get("name", dconfig.get("label"))

            dataset.layout, dataset.trace_params = _dataset_layout(pconfig, dconfig, self.tt_label())
            dataset.dconfig = dconfig

    @staticmethod
    def _set_axis_log_scale(axis):
        axis.type = "log"
        minval = axis.autorangeoptions["minallowed"]
        maxval = axis.autorangeoptions["maxallowed"]
        minval = math.log10(minval) if minval is not None and minval > 0 else None
        maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
        axis.autorangeoptions["minallowed"] = minval
        axis.autorangeoptions["maxallowed"] = maxval

    @staticmethod
    def axis_controlled_by_switches() -> List[str]:
        """
        Return a list of axis names that are controlled by the log10 scale and percentage
        switch buttons, e.g. ["yaxis"]
        """
        return []

    @staticmethod
    def tt_label() -> Optional[str]:
        """Default tooltip label"""
        return None

    def __repr__(self):
        d = {k: v for k, v in self.__dict__.items() if k not in ("datasets", "layout")}
        return f"<{self.__class__.__name__} {self.id} {d}>"

    def add_to_report(self, report) -> str:
        """
        Build and add the plot data to the report, return an HTML wrapper.
        """
        # Setting IDs again now that we have "report" object to guarantee uniqueness
        self.id = report.save_htmlid(self.id)
        for ds in self.datasets:
            ds.uid = self.id
            if len(self.datasets) > 1:  # for flat plots, each dataset will have its own unique ID
                ds.uid = report.save_htmlid(f"{self.id}_{ds.label}", skiplint=True)

        if self.flat:
            html = self.flat_plot()
        else:
            html = self.interactive_plot(report)
            if config.export_plots:
                self.flat_plot()
        return html

    def interactive_plot(self, report) -> str:
        html = '<div class="mqc_hcplot_plotgroup">'

        html += self.__control_panel()

        # This width only affects the space before plot is rendered, and the initial
        # height for the resizing function. For the actual plot container, Plotly will
        # re-calculate the wrapper size after rendering.
        height = self.layout.height
        height_style = f'style="height:{height + 7}px"' if height else ""
        html += f"""
        <div class="hc-plot-wrapper hc-{self.plot_type.value}-wrapper" id="{self.id}-wrapper" {height_style}>
            <div id="{self.id}" class="hc-plot hc-{self.plot_type.value}-plot not_rendered"></div>
            <div class="created-with-multiqc">Created with MultiQC</div>
        </div>"""

        html += "</div>"

        # Saving compressed data for JavaScript to pick up and uncompress.
        dump = self.dump_for_javascript()
        report.plot_data[self.id] = dump
        return html

    def flat_plot(self) -> str:
        html = "".join(
            [
                '<p class="text-info">',
                "<small>" '<span class="glyphicon glyphicon-picture" aria-hidden="true"></span> ',
                "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work ",
                '(see the <a href="https://multiqc.info/docs/development/plots/#interactive--flat-image-plots" target="_blank">docs</a>).',
                "</small>",
                "</p>",
            ]
        )
        html += f'<div class="mqc_mplplot_plotgroup" id="plotgroup-{self.id}" data-pid={self.id}>'

        if not config.simple_output:
            html += self.__control_panel()

        # Go through datasets creating plots
        for ds_idx, dataset in enumerate(self.datasets):
            if self.pconfig.get("save_data_file", True) and self.id != "general_stats_table":
                self.save_data_file(dataset)

            html += self._fig_to_static_html(
                self._make_flat_fig(dataset),
                active=ds_idx == 0 and not self.p_active and not self.l_active,
                uid=dataset.uid if not self.add_log_tab and not self.add_pct_tab else f"{dataset.uid}-cnt",
            )
            if self.add_pct_tab:
                html += self._fig_to_static_html(
                    self._make_flat_fig(dataset, is_pct=True),
                    active=ds_idx == 0 and self.p_active,
                    uid=f"{dataset.uid}-pct",
                )
            if self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_flat_fig(dataset, is_log=True),
                    active=ds_idx == 0 and self.l_active,
                    uid=f"{dataset.uid}-log",
                )
            if self.add_pct_tab and self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_flat_fig(dataset, is_pct=True, is_log=True),
                    active=ds_idx == 0 and self.p_active and self.l_active,
                    uid=f"{dataset.uid}-pct-log",
                )

        html += "</div>"
        return html

    def _btn(self, cls: str, label: str, data_attrs: Dict[str, str] = None, pressed: bool = False) -> str:
        """Build a switch button for the plot."""
        data_attrs = data_attrs.copy() if data_attrs else {}
        if "pid" not in data_attrs:
            data_attrs["pid"] = self.id
        data_attrs = " ".join([f'data-{k}="{v}"' for k, v in data_attrs.items()])
        return f'<button class="btn btn-default btn-sm {cls} {"active" if pressed else ""}" {data_attrs}>{label}</button>\n'

    def buttons(self) -> List[str]:
        """
        Build buttons for control panel
        """
        switch_buttons = ""
        cls = "mpl_switch_group" if self.flat else "interactive-switch-group"
        # Counts / percentages / log10 switches
        if self.add_pct_tab or self.add_log_tab:
            if self.add_pct_tab:
                switch_buttons += self._btn(
                    cls=f"{cls} percent-switch",
                    label=self.pconfig.get("cpswitch_percent_label", "Percentages"),
                    pressed=self.p_active,
                )
            if self.add_log_tab:
                switch_buttons += self._btn(
                    cls=f"{cls} log10-switch",
                    label=self.pconfig.get("logswitch_label", "Log10"),
                    pressed=self.l_active,
                )

        # Buttons to cycle through different datasets
        if len(self.datasets) > 1:
            switch_buttons += f'<div class="btn-group {cls} dataset-switch-group">\n'
            for ds_idx, ds in enumerate(self.datasets):
                data_attrs = {
                    "dataset-index": ds_idx,
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
        if not self.flat:
            export_btn = self._btn(cls="export-plot", label="Export Plot")
        return [switch_buttons, export_btn]

    def __control_panel(self) -> str:
        """
        Add buttons: percentage on/off, log scale on/off, datasets switch panel
        """
        buttons = "\n".join(self.buttons())
        html = f"<div class='row'>\n<div class='col-xs-12'>\n{buttons}\n</div>\n</div>\n\n"
        return html

    def dump_for_javascript(self) -> Dict:
        """Serialise the plot data to pick up in plotly-js"""
        return {
            "id": self.id,
            "layout": self.layout.to_plotly_json(),
            "datasets": [d.dump_for_javascript() for d in self.datasets],
            "plot_type": self.plot_type.value,
            "pct_axis_update": self.pct_axis_update,
            "axis_controlled_by_switches": self._axis_controlled_by_switches,
            "p_active": self.p_active,
            "l_active": self.l_active,
            "square": self.pconfig.get("square"),
            "config": self.pconfig,  # for megaqc
        }

    @abstractmethod
    def save_data_file(self, data: BaseDataset) -> None:
        """
        Save dataset to disk.
        """

    def _make_flat_fig(self, dataset: BaseDataset, is_log=False, is_pct=False) -> go.Figure:
        """
        Create a Plotly Figure object.
        """
        layout = go.Layout(self.layout.to_plotly_json())  # make a copy
        layout.update(**dataset.layout)
        layout.width = layout.width or Plot.FLAT_PLOT_WIDTH
        for axis in self.axis_controlled_by_switches():
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
        return dataset.create_figure(layout, is_log, is_pct)

    @staticmethod
    def _fig_to_static_html(
        fig: go.Figure,
        active: bool,
        uid: str,
    ) -> str:
        """
        Build one static image, return an HTML wrapper.
        """
        assert fig.layout.width
        write_kwargs = dict(
            width=fig.layout.width,  # While interactive plots take full width of screen,
            # for the flat plots we explicitly set width
            height=fig.layout.height,
            scale=2,  # higher detail (retina display)
        )

        # Save the plot to the data directory if export is requested
        if config.export_plots:
            for file_ext in config.export_plot_formats:
                plot_fn = Path(config.plots_dir) / file_ext / f"{uid}.{file_ext}"
                plot_fn.parent.mkdir(parents=True, exist_ok=True)
                if file_ext == "svg":
                    # Cannot add logo to SVGs
                    fig.write_image(plot_fn, **write_kwargs)
                else:
                    img_buffer = io.BytesIO()
                    fig.write_image(img_buffer, **write_kwargs)
                    img_buffer = Plot.add_logo(img_buffer, format=file_ext)
                    with open(plot_fn, "wb") as f:
                        f.write(img_buffer.getvalue())
                    img_buffer.close()

        # Now writing the PNGs for the HTML
        if config.development:
            # Using file written in the config.export_plots block above
            img_src = Path(config.plots_dir_name) / "png" / f"{uid}.png"
        else:
            img_buffer = io.BytesIO()
            fig.write_image(img_buffer, **write_kwargs)
            img_buffer = Plot.add_logo(img_buffer, format="PNG")
            # Convert to a base64 encoded string
            b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
            img_src = f"data:image/png;base64,{b64_img}"
            img_buffer.close()

        # Should this plot be hidden on report load?
        hiding = "" if active else ' style="display:none;"'
        return "".join(
            [
                f'<div class="mqc_mplplot" id="{uid}"{hiding}>',
                f'<img src="{img_src}" height="{fig.layout.height}px" width="{fig.layout.width}px"/>',
                "</div>",
            ]
        )

    @staticmethod
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

    @staticmethod
    def _get_logo() -> str:
        """
        Return the MultiQC logo as a base64 encoded string.
        """
        # logo_path = Path(__file__).parent.parent / "assets" / "img" / "MultiQC_logo.png"
        logo_path = Path(__file__).parent.parent / "assets" / "img" / "MultiQC_logo_dark.svg"
        with open(logo_path, "rb") as f:
            data = f.read()
        return "data:image/svg+xml;base64," + base64.b64encode(data).decode()


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
    pconfig: Dict,
    dconfig: Dict,
    default_tt_label: Optional[str],
) -> Tuple[Dict, Dict]:
    """
    Given plot config and dataset config, set layout and trace params.
    """
    pconfig = pconfig.copy()
    pconfig.update(dconfig)

    # Format on-hover tooltips
    ysuffix = pconfig.get("ysuffix", pconfig.get("tt_suffix"))
    xsuffix = pconfig.get("xsuffix")

    # Options deprecated in 1.21
    ylabformat = pconfig.get("ylab_format", pconfig.get("yLabFormat"))
    if ysuffix is None and ylabformat:
        if "}" in ylabformat:
            ysuffix = ylabformat.split("}")[1]
    xlabformat = pconfig.get("xlab_format", pconfig.get("xLabFormat"))
    if xsuffix is None and xlabformat:
        if "}" in xlabformat:
            xsuffix = xlabformat.split("}")[1]

    # Set or remove space in known suffixes
    KNOWN_SUFFIXES = ["%", "x", "X", "k", "M", " bp", " kbp", " Mbp"]
    for suf in KNOWN_SUFFIXES:
        if ysuffix is not None and ysuffix == suf.strip():
            ysuffix = suf
        if xsuffix is not None and xsuffix == suf.strip():
            xsuffix = suf

    # Set % suffix from ylab if it's in form like "% reads"
    ylab = pconfig.get("ylab")
    xlab = pconfig.get("xlab")
    if ysuffix is None and ylab:
        if "%" in ylab or "percentage" in ylab.lower():
            ysuffix = "%"
        for suf in KNOWN_SUFFIXES:
            if ylab.endswith(f" ({suf.strip()})"):
                ysuffix = suf
    if xsuffix is None and xlab:
        if "%" in xlab or "percentage" in xlab.lower():
            xsuffix = "%"
        for suf in KNOWN_SUFFIXES:
            if xlab.endswith(f" ({suf.strip()})"):
                xsuffix = suf

    if "tt_label" in pconfig:
        # Clean the hover tooltip label, add missing <br> into the beginning, populate suffixes if missing
        tt_label = pconfig["tt_label"]
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
        if ysuffix is not None and "{y}" + ysuffix in tt_label:
            tt_label = tt_label.replace("{y}" + ysuffix, "{y}")
        if xsuffix is not None and "{x}" + xsuffix in tt_label:
            tt_label = tt_label.replace("{x}" + xsuffix, "{x}")

        # add missing line break between the sample name and the key-value pair
        if not tt_label.startswith("<br>"):
            tt_label = "<br>" + tt_label
    else:
        tt_label = default_tt_label

    if tt_label:
        hovertemplate = "<b>%{text}</b>" + tt_label + "<extra></extra>"
    else:
        hovertemplate = None

    # `hoverformat` describes how plain "{y}" or "{x}" are formatted in `hovertemplate`
    y_decimals = pconfig.get("tt_decimals", pconfig.get("y_decimals", pconfig.get("decimalPlaces")))
    y_hoverformat = f",.{y_decimals}f" if y_decimals is not None else None

    x_decimals = pconfig.get("x_decimals")
    x_hoverformat = f",.{x_decimals}f" if x_decimals is not None else None

    layout = dict(
        title=dict(text=pconfig.get("title")),
        xaxis=dict(
            hoverformat=x_hoverformat,
            ticksuffix=xsuffix or "",
            title=dict(text=pconfig.get("xlab")),
            rangemode="tozero" if pconfig.get("xmin") == 0 else "normal",
            autorangeoptions=dict(
                clipmin=pconfig.get("x_clipmin", pconfig.get("xFloor")),
                clipmax=pconfig.get("x_clipmax", pconfig.get("xCeiling")),
                minallowed=pconfig.get("xmin"),
                maxallowed=pconfig.get("xmax"),
            ),
        ),
        yaxis=dict(
            hoverformat=y_hoverformat,
            ticksuffix=ysuffix or "",
            title=dict(text=pconfig.get("ylab")),
            rangemode="tozero" if pconfig.get("ymin") == 0 == 0 else "normal",
            autorangeoptions=dict(
                clipmin=pconfig.get("y_clipmin", pconfig.get("yFloor")),
                clipmax=pconfig.get("y_clipmax", pconfig.get("yCeiling")),
                minallowed=pconfig.get("ymin"),
                maxallowed=pconfig.get("ymax"),
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
