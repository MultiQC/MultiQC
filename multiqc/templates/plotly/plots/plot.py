import base64
import dataclasses
import io
import logging
import random
import string
from abc import abstractmethod, ABC
from enum import Enum
from pathlib import Path
from typing import Dict, Union, List, Optional

import plotly.graph_objects as go

from multiqc.utils import mqc_colour, config

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


@dataclasses.dataclass
class BaseDataset(ABC):
    """
    Structured version of dataset config dictionary
    """

    label: str
    uid: str
    xlab: Optional[str]
    ylab: Optional[str]
    ymax: Optional[float]

    def dump_for_javascript(self) -> Dict:
        """Only dump data added in subclasses. No need to pass all attributes to the dataset"""
        return {k: v for k, v in self.__dict__.items() if k not in BaseDataset.__annotations__.keys()}


class Plot(ABC):
    """Structured version of plot config dictionary"""

    # Default width for flat plots
    FLAT_PLOT_WIDTH = 1100

    def __init__(self, plot_type: PlotType, pconfig: Dict, n_datasets: int):
        if n_datasets == 0:
            raise ValueError("No datasets to plot")

        self.id = pconfig.get("id")
        if self.id is None:  # ID of the plot group
            uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
            is_static_suf = "static_" if config.plots_force_flat else ""
            self.id = f"mqc_{is_static_suf}plot_{uniq_suffix}"

        self.plot_type = plot_type
        self.flat = config.plots_force_flat or (
            not config.plots_force_interactive and n_datasets > config.plots_flat_numseries
        )
        self.pconfig = pconfig
        self.trace_params = dict()

        # Counts / Percentages / Log10 switch
        self.add_log_tab = pconfig.get("logswitch", False) and plot_type in [PlotType.BAR, PlotType.LINE]
        self.add_pct_tab = pconfig.get("cpswitch", True) is not False and plot_type == PlotType.BAR
        self.l_active = self.add_log_tab and pconfig.get("logswitch_active") is True
        self.p_active = self.add_pct_tab and pconfig.get("cpswitch_c_active", True) is not True

        # Per-dataset configurations
        self.datasets: List[BaseDataset] = [
            BaseDataset(
                label=str(i + 1),
                uid=self.id,
                xlab=None,
                ylab=None,
                ymax=None,
            )
            for i in range(n_datasets)
        ]
        data_labels: List[Union[str, Dict[str, str]]] = pconfig.get("data_labels", [])
        for idx, ds in enumerate(self.datasets):
            if idx >= len(data_labels):
                continue
            dl = data_labels[idx]
            if isinstance(dl, str):
                ds.label = dl
            elif isinstance(dl, dict):
                ds.label = dl.get("name", idx + 1)
            else:
                logger.warning(f"Invalid data_labels type: {type(dl)}. Must be a string or a dict.")
            if n_datasets > 1:
                ds.uid += f"_{idx + 1}"
            if isinstance(dl, dict):
                ds.ylab = dl.get("ylab") or ds.label or None
                ds.xlab = dl.get("xlab") or None

        # Format on-hover tooltips
        tt_label = self.tt_label()
        if "tt_label" in self.pconfig:
            tt_label = self.pconfig["tt_label"]
            replace_d = {  # Convert HighCharts format to Plotly format
                "{point.x": "%{x",
                "{point.y": "%{y",
                "{point.category}": "%{x}",
                "<strong>": "<b>",
                "</strong>": "</b>",
                "<br/>": "<br>",
            }
            for k, v in replace_d.items():
                tt_label = tt_label.replace(k, v)
        elif tt_label:
            tt_label += self.pconfig.get("tt_suffix", "")
        if tt_label:
            self.trace_params["hovertemplate"] = "<b>%{text}</b>" + tt_label + "<extra></extra>"

        height = self.pconfig.get("height", 600)
        width = self.pconfig.get("width")
        if self.pconfig.get("square"):
            width = height

        title = self.pconfig.get("table_title", self.pconfig.get("title"))
        if not title:
            logger.error(f"Plot title is not set for {self.id}")

        self.layout = go.Layout(
            title=dict(
                text=title,
                xanchor="center",
                x=0.5,
                font=dict(size=20),
            ),
            xaxis=dict(
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
                title=dict(text=self.pconfig.get("xlab") or (self.datasets[0].xlab if self.datasets else None)),
                rangemode="tozero" if self.pconfig.get("xmin") == 0 else "normal",
                range=[self.pconfig.get("xmin"), self.pconfig.get("xmax")],
                tickformat=".0%" if self.p_active and self.add_pct_tab else None,
            ),
            yaxis=dict(
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
                title=dict(text=self.pconfig.get("ylab") or (self.datasets[0].ylab if self.datasets else None)),
                rangemode="tozero" if self.pconfig.get("ymin") == 0 else "normal",
                range=[self.pconfig.get("ymin"), self.pconfig.get("ymax", self.pconfig.get("yCeiling"))],
                # Default precision for floating numbers is too high - allowing to override it
                hoverformat=f".{pconfig['tt_decimals']}f" if "tt_decimals" in pconfig else None,
            ),
            height=height,
            width=width,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(color="Black", family="Lucida Grande"),
            colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
            autosize=True,
            margin=dict(
                pad=10,  # pad sample names in a bar graph a bit
                t=50,  # more compact title
                r=10,  # remove excessive whitespace on the right
            ),
            hoverlabel=dict(
                namelength=-1,  # do not crop sample names inside hover label <extra></extra>
            ),
            modebar=dict(
                bgcolor="rgba(0, 0, 0, 0)",
                color="rgba(0, 0, 0, 0.5)",
                activecolor="rgba(0, 0, 0, 1)",
            ),
            showlegend=self.flat,
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.15,
                xanchor="center",
                x=0.5,
            )
            if self.flat
            else None,
        )
        if self.add_log_tab and self.l_active:
            for axis in self.axis_controlled_by_switches():
                self.layout[axis].type = "log"

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
        return html

    def interactive_plot(self, report) -> str:
        html = '<div class="mqc_hcplot_plotgroup">'

        html += self._control_panel()

        height_style = f'style="height:{self.layout.height}px"' if self.layout.height else ""
        html += f"""
        <div class="hc-plot-wrapper hc-{self.plot_type.value}-wrapper" id="{self.id}-wrapper" {height_style}>
            <div id="{self.id}" class="hc-plot hc-{self.plot_type.value}-plot not_rendered"></div>
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
                '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).',
                "</small>",
                "</p>",
            ]
        )
        html += f'<div class="mqc_mplplot_plotgroup" id="plotgroup-{self.id}" data-pid={self.id}>'

        if not config.simple_output:
            html += self._control_panel()

        # Go through datasets creating plots
        for ds_idx, dataset in enumerate(self.datasets):
            if self.pconfig.get("save_data_file", True):
                self.save_data_file(dataset)

            html += self._fig_to_static_html(
                self._make_fig(dataset),
                active=ds_idx == 0 and not self.p_active and not self.l_active,
                uid=dataset.uid if not self.add_log_tab and not self.add_pct_tab else f"{dataset.uid}-cnt",
            )
            if self.add_pct_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_pct=True),
                    active=ds_idx == 0 and self.p_active,
                    uid=f"{dataset.uid}-pct",
                )
            if self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_log=True),
                    active=ds_idx == 0 and self.l_active,
                    uid=f"{dataset.uid}-log",
                )
            if self.add_pct_tab and self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_pct=True, is_log=True),
                    active=ds_idx == 0 and self.p_active and self.l_active,
                    uid=f"{dataset.uid}-pct-log",
                )

        html += "</div>"
        return html

    def _btn(self, cls: str, label: str, data_attrs: Dict[str, str] = None, pressed: bool = False) -> str:
        """Build a switch button for the plot."""
        data_attrs = data_attrs.copy() if data_attrs else {}
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
                    "ylab": ds.ylab,
                    "ymax": ds.ymax,
                    "xlab": ds.xlab,
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

    def _control_panel(self) -> str:
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
            "plot_type": self.plot_type.value,
            "layout": self.layout.to_plotly_json(),
            "trace_params": self.trace_params,
            "axis_controlled_by_switches": self.axis_controlled_by_switches(),
            "datasets": [d.dump_for_javascript() for d in self.datasets],
            "p_active": self.p_active,
            "l_active": self.l_active,
            "square": self.pconfig.get("square"),
        }

    @abstractmethod
    def create_figure(self, layout: go.Layout, dataset: BaseDataset, is_log=False, is_pct=False):
        """
        To be overridden by specific plots: create a Plotly figure for a dataset, update layout if needed.
        """
        pass

    @abstractmethod
    def save_data_file(self, data: BaseDataset) -> None:
        """
        Save dataset to disk.
        """

    def _make_fig(self, dataset: BaseDataset, is_log=False, is_pct=False) -> go.Figure:
        """
        Create a Plotly Figure object.
        """
        layout = go.Layout(**self.layout.to_plotly_json())  # make a copy
        layout.width = layout.width or Plot.FLAT_PLOT_WIDTH
        for axis in self.axis_controlled_by_switches():
            layout[axis].type = "log" if is_log else "linear"
            if is_pct:
                layout[axis].tickformat = ".0%"
        return self.create_figure(layout, dataset, is_log, is_pct)

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
                fig.write_image(plot_fn, **write_kwargs)

        # Now writing the PNGs for the HTML
        write_kwargs["format"] = "png"

        if config.development:
            img_src = Path(config.plots_dir_name) / "png" / f"{uid}.png"
        else:
            # Output the figure to a base64 encoded string
            img_buffer = io.BytesIO()
            fig.write_image(img_buffer, **write_kwargs)
            b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
            img_buffer.close()
            img_src = f"data:image/png;base64,{b64_img}"

        # Should this plot be hidden on report load?
        hiding = "" if active else ' style="display:none;"'
        return "".join(
            [
                f'<div class="mqc_mplplot" id="{uid}"{hiding}>',
                f'<img src="{img_src}" height="{fig.layout.height}px" width="{fig.layout.width}px"/>',
                "</div>",
            ]
        )
