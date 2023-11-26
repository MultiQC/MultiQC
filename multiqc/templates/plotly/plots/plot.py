# from abc import ABC, abstractmethod
# from typing import Dict, List
#
# import plotly.graph_objects as go
#
# from multiqc.utils import config, report
import logging
from collections import namedtuple
from pprint import pprint
from typing import Dict, Union, List, Optional
import plotly.graph_objects as go

from multiqc.utils import mqc_colour, config

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

logger = logging.getLogger(__name__)


class PlotSettings:
    """
    Build a structured object from the plot configuration dictionary.
    """

    def __init__(self, pconfig: Dict, num_datasets: int):
        self.title = pconfig.get("title")
        self.id = pconfig.get("id")
        self.height = pconfig.get("height")
        self.save_data_file: bool = pconfig.get("save_data_file", True)

        # Counts / Percentages / Log10 switch
        self.add_pct_tab = pconfig.get("cpswitch", True) is not False
        self.add_log_tab = pconfig.get("logswitch", False)
        self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.l_label = pconfig.get("logswitch_label", "Log10")
        self.c_active = "active"
        self.p_active = ""
        self.l_active = ""
        self.stacking = pconfig.get("stacking", "stack")
        if self.add_pct_tab and not config.simple_output:
            if pconfig.get("logswitch_active") is True:
                self.l_active = "active"
            elif pconfig.get("cpswitch_c_active", True) is True:
                self.c_active = "active"
            else:
                self.p_active = "active"
                self.stacking = "relative"

        self.categories: Optional[List[str]] = pconfig.get("categories")

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

        self.xmax = pconfig.get("xmax")
        self.ymax = pconfig.get("ymax")
        self.xmin = pconfig.get("xmin")
        self.ymin = pconfig.get("ymin")

    def __repr__(self):
        return f"<{self.__class__.__name__} {pprint(self.__dict__)}>"


def base_layout(settings) -> go.Layout:
    """
    Layout object for the line plot.
    """
    return go.Layout(
        title=dict(
            text=settings.title,
            xanchor="center",
            x=0.5,
            font=dict(size=20),
        ),
        yaxis=dict(
            title=dict(text=settings.ylab),
            gridcolor="rgba(0,0,0,0.1)",
        ),
        xaxis=dict(
            title=dict(text=settings.xlab),
            gridcolor="rgba(0,0,0,0.1)",
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
