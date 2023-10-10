from abc import ABC, abstractmethod
from typing import Dict, List

import plotly.graph_objects as go

from multiqc.utils import config, report


class AbstractPlot(ABC):
    def __init__(self, pconfig):
        self.title = pconfig.get("title")
        self._id = pconfig.get("id")
        self.xlab = pconfig.get("xlab")
        self.ylab = pconfig.get("ylab")
        self.height = pconfig.get("height")

        self.add_pct_tab = pconfig.get("cpswitch", True) is not False
        self.add_log_tab = pconfig.get("logswitch", False)

        # Counts / Percentages Switch
        self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.l_label = pconfig.get("logswitch_label", "Log10")
        self.c_active = ""
        self.p_active = ""
        self.l_active = ""
        self.stacking = pconfig.get("stacking", "normal")
        if self.add_pct_tab and not config.simple_output:
            if pconfig.get("logswitch_active") is True:
                self.l_active = "active"
            elif pconfig.get("cpswitch_c_active", True) is True:
                self.c_active = "active"
            else:
                self.p_active = "active"
                self.stacking = "percent"

        self.data_labels: List = pconfig.get("data_labels", [])
        self.save_data_file = pconfig.get("save_data_file", True)

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, id: str):
        # Sanitise plot ID and check for duplicates
        self._id = report.save_htmlid(id)

    @abstractmethod
    def _layout(self) -> go.Layout:
        """
        Make a Plotly Layout object. Will be serialised and used by plotly.js to create
        an interactive plot, or used to create a static version if flat plots are
        requested.
        """
        pass

    @abstractmethod
    def _figure(
        self,
        layout: go.Layout,
        data_by_cat: List[Dict],
        sample_names: List[str],
    ) -> go.Figure:
        """
        Make a Plotly Figure object. Will be called multiple times if a multi-tab plot
        is requested. Used to build flat versions. For an interactive plot, a similar
        function will be called in JavaScript based on plotly.js.
        """
        pass

    @abstractmethod
    def add_to_report(
        self,
        data_by_cat_lists: List[List[Dict]],
        samples_lists: List[List[str]],
    ) -> str:
        """
        Build and add the plot data to the report, return an HTML wrapper.
        :param data_by_cat_lists: List of lists of dicts with the keys:
            {name, color, data}, where `name` is the category name, `color` is the
            color of the bar, and `data` is a list of values for each sample.
            Each outer list will correspond a separate tab.
        :param samples_lists: List of lists of bar names (i.e., sample names).
            Similarly, each outer list will correspond to a separate tab.
        :return: an HTML string
        """
        pass

    @staticmethod
    def _add_to_report_interactive(
        self,
        data_by_cat_lists: List[List[Dict]],
        samples_lists: List[List[str]],
    ) -> str:
        """
        Add data to the report and build an HTML wrapper for interactive plotly-js.
        """
        pass

    @abstractmethod
    def _add_to_report_flat(
        self,
        data_by_cat_lists: List[List[Dict]],
        samples_lists: List[List[str]],
    ) -> str:
        """
        Build HTML with plots as static images.
        """
        pass
