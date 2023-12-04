# from abc import ABC, abstractmethod
# from typing import Dict, List
#
# import plotly.graph_objects as go
#
# from multiqc.utils import config, report
import base64
import io
import logging
import random
import string
from abc import abstractmethod, ABC
from enum import Enum
from pathlib import Path
from pprint import pprint
from typing import Dict, Union, List

import plotly.graph_objects as go

from multiqc.plots.bargraph import get_template_mod
from multiqc.utils import mqc_colour, config

logger = logging.getLogger(__name__)


class PlotType(Enum):
    LINE = "xy_line"
    BAR = "bar_graph"
    SCATTER = "scatter"
    HEATMAP = "heatmap"


class Plot(ABC):
    """Structured version of plot config dictionary"""

    def __init__(self, plot_type: PlotType, pconfig: Dict, num_datasets: int):
        self.id = pconfig.get("id")
        if self.id is None:  # ID of the plot group
            uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
            is_static_suf = "static_" if config.plots_force_flat else ""
            self.id = f"mqc_{is_static_suf}plot_{uniq_suffix}"

        self.plot_type = plot_type
        self.title = pconfig.get("title")
        self.height = pconfig.get("height")
        self.do_save_data_file: bool = pconfig.get("save_data_file", True)

        # Counts / Percentages / Log10 switch
        self.add_log_tab = pconfig.get("logswitch", False) and plot_type in [PlotType.BAR, PlotType.LINE]
        self.add_pct_tab = pconfig.get("cpswitch", True) is not False and plot_type == PlotType.BAR
        self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
        self.l_label = pconfig.get("logswitch_label", "Log10")
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.l_active = self.add_log_tab and pconfig.get("logswitch_active") is True
        self.p_active = self.add_pct_tab and pconfig.get("cpswitch_c_active", True) is not True

        # Per-dataset configurations
        data_labels: List[Union[str, Dict[str, str]]] = pconfig.get("data_labels", [])
        for idx in range(num_datasets):
            if len(data_labels) <= idx:
                data_labels.append({})
            dl = data_labels[idx]
            if isinstance(dl, str):
                data_labels[idx] = {"name": dl}
            elif isinstance(dl, dict):
                data_labels[idx]["name"] = dl.get("name", idx + 1)
                data_labels[idx]["ylab"] = dl.get("ylab") or dl.get("name") or None
            else:
                logger.warning(
                    f"Invalid data_labels type: {type(data_labels[idx])}. " f"Must be a string or a dictionary."
                )
                data_labels[idx] = {"name": str(idx + 1)}
            data_labels[idx]["uid"] = self.id
            if num_datasets > 1:
                data_labels[idx]["uid"] += f"_{idx + 1}"
        self.data_labels: List[Dict[str, str]] = data_labels

        # Add initial axis labels if defined in `data_labels` but not main config
        self.ylab = pconfig.get("ylab") or self.data_labels[0].get("ylab")
        self.xlab = pconfig.get("xlab") or self.data_labels[0].get("xlab")
        self.xmin = pconfig.get("xmin")
        self.xmax = pconfig.get("xmax")
        self.ymin = pconfig.get("ymin")
        self.ymax = pconfig.get("ymax")

        # Default state to dump to JSON to load with plotly.js
        self.config = {
            "p_active": self.p_active,
            "l_active": self.l_active,
        }

    def __repr__(self):
        return f"<{self.__class__.__name__} {pprint(self.__dict__)}>"

    def layout(self) -> go.Layout:
        """
        Layout object for the line plot.
        """
        layout = go.Layout(
            title=dict(
                text=self.title,
                xanchor="center",
                x=0.5,
                font=dict(size=20),
            ),
            xaxis=dict(
                title=dict(text=self.xlab),
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
                rangemode="tozero" if self.xmin == 0 else "normal",
                range=[self.xmin, self.xmax],
            ),
            yaxis=dict(
                title=dict(text=self.ylab),
                gridcolor="rgba(0,0,0,0.1)",
                zerolinecolor="rgba(0,0,0,0.1)",
                rangemode="tozero" if self.ymin == 0 else "normal",
                range=[self.ymin, self.ymax],
            ),
            height=self.height,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            font={"color": "Black", "family": "Lucida Grande"},
            colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
            showlegend=False,
            autosize=True,
            margin_pad=10,  # pad sample names a bit
            hoverlabel_namelength=-1,  # do not crop sample names in hover labels
            annotations=[
                dict(
                    text="Created with MultiQC",
                    font=dict(size=10, color="rgba(0,0,0,0.5)"),
                    xanchor="right",
                    yanchor="top",
                    x=1.05,
                    y=1.05,
                    textangle=-90,
                    # yanchor="bottom",
                    # x=1.07,
                    # y=-0.35,
                    xref="paper",
                    yref="paper",
                    showarrow=False,
                )
            ],
            modebar=dict(
                bgcolor="rgba(0, 0, 0, 0)",
                color="rgba(0, 0, 0, 0.5)",
                activecolor="rgba(0, 0, 0, 1)",
            ),
        )
        return layout

    def add_to_report(self, datasets: List[List[Dict]], report) -> str:
        """
        Build and add the plot data to the report, return an HTML wrapper.
        """
        self.id = report.save_htmlid(self.id)
        for dl in self.data_labels:
            dl["uid"] = self.id
            if len(datasets) > 1:  # for flat plots, each dataset will have its own unique ID
                dl["uid"] = report.save_htmlid(f"{self.id}_{dl['name']}", skiplint=True)

        if config.plots_force_flat:
            html = self.flat_plot(datasets)
        else:
            html = self.interactive_plot(report, datasets)
        return html

    def flat_plot(self, datasets: List[List]) -> str:
        html = (
            '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
            + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
            + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
        )
        html += f'<div class="mqc_mplplot_plotgroup" id="{self.id}">'

        if not config.simple_output:
            html += self._buttons(datasets, cls="mpl_switch_group")

        # Go through datasets creating plots
        for ds_idx, dataset in enumerate(datasets):
            uid = self.data_labels[ds_idx]["uid"]

            if self.do_save_data_file:
                self.save_data_file(dataset, uid=uid)

            html += self._fig_to_static_html(
                self._make_fig(dataset),
                active=ds_idx == 0 and not self.p_active and not self.l_active,
                uid=uid,
            )
            if self.add_pct_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_pct=True),
                    active=ds_idx == 0 and self.p_active,
                    uid=f"{uid}_pc",
                )
            if self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_log=True),
                    active=ds_idx == 0 and self.l_active,
                    uid=f"{uid}_log",
                )
            if self.add_pct_tab and self.add_log_tab:
                html += self._fig_to_static_html(
                    self._make_fig(dataset, is_pct=True, is_log=True),
                    active=ds_idx == 0 and self.p_active and self.l_active,
                    uid=f"{uid}_pc_log",
                )

        html += "</div>"
        return html

    def _make_fig(self, dataset: List, is_log=False, is_pct=False) -> go.Figure:
        """
        Create a Plotly Figure object.
        """
        layout = self.layout()
        if is_pct:
            layout.update(
                {
                    "xaxis.tickformat": ".0%",
                    # "xaxis.title.text": self.p_label,
                }
            )
        if is_log:
            layout.update(
                {
                    "xaxis.type": "log",
                    # "xaxis.title.text": layout.xaxis.title.text + " " + self.l_label
                }
            )
        fig = go.Figure(layout=layout)
        self.populate_figure(fig, dataset, is_log, is_pct)
        return fig

    @abstractmethod
    def populate_figure(self, fig: go.Figure, dataset: List, is_log=False, is_pct=False):
        """
        To be overridden by specific plots: add traces, updated layout if needed.
        """
        pass

    @abstractmethod
    def save_data_file(self, dataset: List, uid: str) -> None:
        """
        Save dataset to disk.
        """

    def _buttons(self, datasets, cls=""):
        """
        Add buttons: percentage on/off, log scale on/off, datasets switch panel
        """

        def _btn(cls: str, pid: str, active: bool, label: str, data_attrs: Dict[str, str] = None) -> str:
            """Build a switch button for the plot."""
            data_attrs = data_attrs.copy() if data_attrs else {}
            data_attrs["pid"] = pid
            data_attrs = " ".join([f'data-{k}="{v}"' for k, v in data_attrs.items()])
            return f'<button class="btn btn-default btn-sm {cls} {"active" if active else ""}" {data_attrs}>{label}</button>\n'

        html = ""
        # Counts / percentages / log10 switches
        if self.add_pct_tab or self.add_log_tab:
            if self.add_pct_tab:
                html += _btn(
                    cls=f"{cls} percent-switch",
                    pid=self.id,
                    active=self.p_active,
                    label=self.p_label,
                )
            if self.add_log_tab:
                html += _btn(
                    cls=f"{cls} log10-switch",
                    pid=self.id,
                    active=self.l_active,
                    label=self.l_label,
                )
            if len(datasets) > 1:
                html += " &nbsp; &nbsp; "

        # Buttons to cycle through different datasets
        if len(datasets) > 1:
            html += f'<div class="btn-group {cls} dataset-switch-group">\n'
            for ds_idx, ds in enumerate(datasets):
                dl: Dict[str, str] = self.data_labels[ds_idx]
                data_attrs = {k: dl[k] for k in ["ylab", "ymax", "xlab"] if k in dl}
                data_attrs["dataset-index"] = ds_idx
                # For flat plots, we will generate separate flat images for each
                # dataset and view, so have to save individual image IDs.
                data_attrs["dataset-uid"] = dl["uid"]
                html += _btn(
                    cls="",
                    pid=self.id,
                    active=ds_idx == 0,
                    label=dl["name"],
                    data_attrs=data_attrs,
                )
            html += "</div>\n\n"
        return html

    def serialise(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        return {
            "id": self.id,
            "plot_type": self.plot_type.value,
            "layout": self.layout().to_plotly_json(),
            # TODO: save figures to JSON, not datasets?
            # "figures": [self._make_fig(dataset).to_plotly_json() for dataset in datasets],
            "config": self.config,
        }

    def interactive_plot(self, report, datasets) -> str:
        html = '<div class="mqc_hcplot_plotgroup">'
        html += self._buttons(datasets, cls="interactive-switch-group")

        # Plot HTML
        html += """
        <div class="hc-plot-wrapper"{height}>
            <div id="{id}" class="hc-plot not_rendered hc-{plot_type}-plot"></div>
        </div>""".format(
            id=self.id,
            height=f' style="height:{self.height}px"' if self.height else "",
            plot_type=self.plot_type,
        )
        html += "</div>"

        # Saving compressed data for JavaScript to pick up and uncompress.
        d = self.serialise()
        d["datasets"] = datasets
        report.plot_data[self.id] = d
        return html

    @staticmethod
    def _fig_to_static_html(
        fig: go.Figure,
        active: bool,
        uid: str,
    ) -> str:
        """
        Build one static image, return an HTML wrapper.
        """
        # Save the plot to the data directory if export is requested
        if config.export_plots:
            for file_ext in config.export_plot_formats:
                # Make the directory if it doesn't already exist
                plot_dir = Path(config.plots_dir) / file_ext
                if not plot_dir.exists():
                    plot_dir.mkdir(parents=True, exist_ok=True)
                # Save the plot
                plot_fn = Path(plot_dir) / f"{uid}.{file_ext}"
                fig.write_image(
                    plot_fn,
                    format=file_ext,
                    width=fig.layout.width,
                    height=fig.layout.height,
                    scale=1,
                )

        fig_write_args = dict(
            format="png",
            width=1100,
            height=fig.layout.height,
            scale=1,
        )

        # Output the figure to a base64 encoded string
        if getattr(get_template_mod(), "base64_plots", True) is True:
            img_buffer = io.BytesIO()
            fig.write_image(img_buffer, **fig_write_args)
            b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
            img_buffer.close()
            img_src = f"data:image/png;base64,{b64_img}"

        # Link to the saved image
        else:
            plot_relpath = Path(config.plots_dir_name) / "png" / f"{uid}.png"
            plot_relpath.parent.mkdir(parents=True, exist_ok=True)
            fig.write_image(plot_relpath, **fig_write_args)
            img_src = plot_relpath

        # Should this plot be hidden on report load?
        hide_div = "" if active else ' style="display:none;"'
        return f'<div class="mqc_mplplot" id="flat-{uid}"{hide_div}><img src="{img_src}" /></div>'
