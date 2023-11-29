# from abc import ABC, abstractmethod
# from typing import Dict, List
#
# import plotly.graph_objects as go
#
# from multiqc.utils import config, report
import logging
import random
import string
from abc import abstractmethod, ABC
from collections import namedtuple
from pprint import pprint
from typing import Dict, Union, List
import plotly.graph_objects as go

from multiqc.utils import mqc_colour, config

logger = logging.getLogger(__name__)


View = namedtuple(
    "View",
    [
        "data",
        "active",
        "suffix",
        "label",
        "xaxis_tickformat",
    ],
)


class Plot(ABC):
    def __init__(self, plot_type: str, pconfig: Dict, num_datasets: int):
        self.plot_type = plot_type
        self.title = pconfig.get("title")
        self.id = pconfig.get("id")
        self.height = pconfig.get("height")
        self.save_data_file: bool = pconfig.get("save_data_file", True)

        # Counts / Percentages / Log10 switch
        self.add_log_tab = pconfig.get("logswitch", False)
        self.add_pct_tab = False
        self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
        self.l_label = pconfig.get("logswitch_label", "Log10")
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.c_active = "active"
        self.l_active = ""
        self.p_active = ""
        if self.add_log_tab and not config.simple_output and pconfig.get("logswitch_active") is True:
            self.c_active = ""
            self.l_active = "active"

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
        self.data_labels: List[Dict[str, str]] = data_labels

        # Add initial axis labels if defined in `data_labels` but not main config
        self.ylab = pconfig.get("ylab") or self.data_labels[0].get("ylab")
        self.xlab = pconfig.get("xlab") or self.data_labels[0].get("xlab")
        self.xmin = pconfig.get("xmin")
        self.xmax = pconfig.get("xmax")
        self.ymin = pconfig.get("ymin")
        self.ymax = pconfig.get("ymax")

    def __repr__(self):
        return f"<{self.__class__.__name__} {pprint(self.__dict__)}>"

    def layout(self) -> go.Layout:
        """
        Layout object for the line plot.
        """
        # if pconfig.xmin is not None or pconfig.xmax is not None:
        #     xrange = [pconfig.xmin, pconfig.xmax]
        # if pconfig.ymin is not None or pconfig.ymax is not None:
        #     yrange = [pconfig.ymin, pconfig.ymax]
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
                rangemode="tozero" if self.xmin == 0 else "normal",
                range=[self.xmin, self.xmax],
            ),
            yaxis=dict(
                title=dict(text=self.ylab),
                gridcolor="rgba(0,0,0,0.1)",
                rangemode="tozero" if self.ymin == 0 else "normal",
                range=[self.ymin, self.ymax],
            ),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            font={"color": "Black", "family": "Lucida Grande"},
            colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
            showlegend=False,
            autosize=True,
            margin=dict(
                pad=10,  # pad sample names a bit
            ),
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

    def add_to_report(
        self,
        datasets: List[List[Dict]],
        report,
    ) -> str:
        """
        Build and add the plot data to the report, return an HTML wrapper.
        """
        uniq_suffix = "".join(random.sample(string.ascii_lowercase, 10))
        is_static_suf = "static_" if config.plots_force_flat else ""
        pid = report.save_htmlid(self.id)
        if pid is None:  # ID of the plot group
            pid = report.save_htmlid(f"mqc_{is_static_suf}plot_{uniq_suffix}")

        if config.plots_force_flat:
            html = self.flat_plot(
                datasets,
                pid,
                # Can't use interactivity, so we will have to generate separate flat images for each
                # dataset and view. So have to make sure the individual image IDs are unique across the report:
                pids=[report.save_htmlid(f"{pid}_{dl['name']}", skiplint=True) for dl in self.data_labels],
            )
        else:
            html = self.interactive_plot(
                report=report,
                datasets=datasets,
                pid=pid,
                layout=self.layout(),
            )
        return html

    def flat_plot(
        self,
        datasets: List[List],
        pid: str,
        pids: List[str],
    ) -> str:
        html = (
            '<p class="text-info"><small><span class="glyphicon glyphicon-picture" aria-hidden="true"></span> '
            + "Flat image plot. Toolbox functions such as highlighting / hiding samples will not work "
            + '(see the <a href="http://multiqc.info/docs/#flat--interactive-plots" target="_blank">docs</a>).</small></p>'
        )
        html += f'<div class="mqc_mplplot_plotgroup" id="{pid}">'

        # Counts / Percentages Switch
        if (self.add_pct_tab or self.add_log_tab) and not config.simple_output:
            html += (
                '<div class="btn-group mpl_switch_group flat_plot_toggle_percent_log"> \n'
                + f'<button class="btn btn-default btn-sm {self.c_active} counts">{self.c_label}</button> \n'
            )
            if self.add_pct_tab:
                html += f'<button class="btn btn-default btn-sm {self.p_active} percents">{self.p_label}</button> \n'
            if self.add_log_tab:
                html += f'<button class="btn btn-default btn-sm {self.l_active} log10">{self.l_label}</button> \n'
            # if self.add_pct_tab and self.add_log_tab:
            #     html += f'<button class="btn btn-default btn-sm {self.l_active} {self.p_active} log10_percents">{self.l_label} {self.p_label}</button> \n'
            html += "</div>"
            if len(datasets) > 1:
                html += " &nbsp; &nbsp; "

        # Buttons to cycle through different datasets
        if len(datasets) > 1 and not config.simple_output:
            html += '<div class="btn-group mpl_switch_group flat_plot_switch_datasets">\n'
            for pidx, ds in enumerate(datasets):
                pid = pids[pidx]
                active = "active" if pidx == 0 else ""
                name = self.data_labels[pidx]["name"]
                html += f'<button class="btn btn-default btn-sm {active}" data-target="#{pid}">{name}</button>\n'
            html += "</div>\n\n"

        # Go through datasets creating plots
        for pidx, (pid, dataset) in enumerate(zip(pids, datasets)):
            html += self._flat_imgs_for_dataset(pidx, pid, dataset)

        html += "</div>"
        return html

    @abstractmethod
    def _flat_imgs_for_dataset(
        self,
        pidx: int,
        pid: str,
        dataset: List,
    ) -> str:
        pass

    def interactive_plot(
        self,
        report,
        datasets,
        pid: str,
        layout: go.Layout,
    ) -> str:
        def switch_button(cls: str, pid: str, active: str, label: str, data_attrs: str = "") -> str:
            """Build a switch button for the plot."""
            return (
                f'<button class="{cls} btn btn-default btn-sm {active}" '
                f'data-target="{pid}" '
                f"{data_attrs}>"
                f"{label}"
                "</button>\n"
            )

        html = '<div class="mqc_hcplot_plotgroup">'
        # Counts / Percentages / Log Switches
        if self.add_pct_tab or self.add_log_tab:
            if self.add_pct_tab:
                html += switch_button(
                    "switch_percent",
                    pid,
                    self.p_active,
                    self.p_label,
                )
            if self.add_log_tab:
                html += switch_button(
                    "switch_log10",
                    pid,
                    self.l_active,
                    self.l_label,
                )
            if len(datasets) > 1:
                html += " &nbsp; &nbsp; "

        # Buttons to cycle through different datasets
        if len(datasets) > 1:
            html += '<div class="btn-group dataset_switch_group">\n'
            for k, ds in enumerate(datasets):
                active = "active" if k == 0 else ""
                dl: Dict[str, str] = self.data_labels[k]
                name = dl["name"]
                ylab = f'data-ylab="{dl["ylab"]}"' if "ylab" in dl else ""
                ymax = f'data-ylab="{dl["ymax"]}"' if "ymax" in dl else ""
                xlab = f'data-ylab="{dl["xlab"]}"' if "xlab" in dl else ""
                html += switch_button(
                    "",
                    pid,
                    active,
                    name,
                    f'{ylab} {ymax} {xlab} data-dataset_index="{k}"',
                )
            html += "</div>\n\n"

        # Plot HTML
        html += """
        <div class="hc-plot-wrapper"{height}>
            <div id="{id}" class="hc-plot not_rendered hc-{plot_type}-plot"></div>
        </div>""".format(
            id=pid,
            height=f' style="height:{self.height}px"' if self.height else "",
            plot_type=self.plot_type,
        )
        # Close wrapping div
        html += "</div>"

        # Saving compressed data for JavaScript to pick up and uncompress.
        report.num_hc_plots += 1
        report.plot_data[pid] = {
            "plot_type": self.plot_type,
            "datasets": datasets,
            "layout": layout.to_plotly_json(),
            "pconfig": self.__dict__,
        }
        # Parameters to be toggled by switch buttons:
        report.plot_data[pid]["active_dataset_idx"] = 0
        if self.add_log_tab:
            report.plot_data[pid]["l_active"] = self.l_active
        if self.add_pct_tab:
            report.plot_data[pid]["p_active"] = self.p_active
        return html


# class AbstractPlot(ABC):
#     def __init__(self, pconfig):
#         self.title = pconfig.get("title")
#         self._id = pconfig.get("id")
#         self.xlab = pconfig.get("xlab")
#         self.ylab = pconfig.get("ylab")
#         self.height = pconfig.get("height")
#
#         self.add_pct_tab = pconfig.get("cpswitch", True) is not False
#         self.add_log_tab = pconfig.get("logswitch", False)
#
#         # Counts / Percentages Switch
#         self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
#         self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
#         self.l_label = pconfig.get("logswitch_label", "Log10")
#         self.c_active = ""
#         self.p_active = ""
#         self.l_active = ""
#         self.stacking = pconfig.get("stacking", "normal")
#         if self.add_pct_tab and not config.simple_output:
#             if pconfig.get("logswitch_active") is True:
#                 self.l_active = "active"
#             elif pconfig.get("cpswitch_c_active", True) is True:
#                 self.c_active = "active"
#             else:
#                 self.p_active = "active"
#                 self.stacking = "percent"
#
#         self.data_labels: List = pconfig.get("data_labels", [])
#         self.save_data_file = pconfig.get("save_data_file", True)
#
#     @property
#     def id(self) -> str:
#         return self._id
#
#     @id.setter
#     def id(self, id: str):
#         # Sanitise plot ID and check for duplicates
#         self._id = report.save_htmlid(id)
#
#     @abstractmethod
#     def _layout(self) -> go.Layout:
#         """
#         Make a Plotly Layout object. Will be serialised and used by plotly.js to create
#         an interactive plot, or used to create a static version if flat plots are
#         requested.
#         """
#         pass
#
#     @abstractmethod
#     def _figure(self, layout: go.Layout, *args) -> go.Figure:
#         """
#         Make a Plotly Figure object. Will be called multiple times if a multi-tab plot
#         is requested. Used to build flat versions. For an interactive plot, a similar
#         function will be called in JavaScript based on plotly.js.
#         """
#         pass
#
#     @abstractmethod
#     def add_to_report(self, *args) -> str:
#         """
#         Build and add the plot data to the report, return an HTML wrapper.
#         :param data_by_cat_lists: List of lists of dicts with the keys:
#             {name, color, data}, where `name` is the category name, `color` is the
#             color of the bar, and `data` is a list of values for each sample.
#             Each outer list will correspond a separate tab.
#         :param samples_lists: List of lists of bar names (i.e., sample names).
#             Similarly, each outer list will correspond to a separate tab.
#         :return: an HTML string
#         """
#         pass
#
#     @staticmethod
#     def _add_to_report_interactive(self, *args) -> str:
#         """
#         Add data to the report and build an HTML wrapper for interactive plotly-js.
#         """
#         pass
#
#     @abstractmethod
#     def _add_to_report_flat(self, *args) -> str:
#         """
#         Build HTML with plots as static images.
#         """
#         pass
