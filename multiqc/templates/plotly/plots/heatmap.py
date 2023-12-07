import logging
from typing import Dict, List, Union
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, Dataset

logger = logging.getLogger(__name__)


ElemT = Union[str, float, int]


def plot(
    dataset: List[List[ElemT]],
    xcats: List[str],
    ycats: List[str],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param dataset: One dataset. A dataset is a list of rows of values
    :param xcats: Labels for X axis
    :param ycats: Labels for Y axis
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = HeatmapPlot(pconfig, [dataset], xcats, ycats)

    from multiqc.utils import report

    return p.add_to_report(report)


class HeatmapPlot(Plot):
    def __init__(self, pconfig: Dict, datasets: List, xcats: List[str], ycats: List[str]):
        super().__init__(PlotType.HEATMAP, pconfig, datasets)
        self.xcats = xcats
        self.ycats = ycats

        # Determining the size of the plot to reasonably display data
        # without cluttering too much
        n = len(self.xcats)
        MIN_PLOT_HEIGHT = 500
        MAX_PLOT_HEIGHT = 2560
        height = MIN_PLOT_HEIGHT
        font_size = 12
        nticks = n  # Making sure all ticks are displayed
        if n < 15:
            pass
        elif n < 20:
            height = max(height, n * 30)
        elif n < 30:
            height = max(height, n * 26)
        elif n < 50:
            height = max(height, n * 18)
            font_size = 10
        else:
            height = max(height, n * 13)
            font_size = 8

        if height > MAX_PLOT_HEIGHT:
            # Cap and allow skipping ticks at this point
            height = MAX_PLOT_HEIGHT
            nticks = None

        self.layout.xaxis.tickangle = 45
        self.layout.font.size = font_size
        self.layout.xaxis.nticks = self.layout.yaxis.nticks = nticks
        self.layout.height = self.layout.height or height
        self.layout.width = self.layout.width or height

    def serialise(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().serialise()
        d["xcats"] = self.xcats
        d["ycats"] = self.ycats
        d["xcats_samples"] = self.pconfig.get("xcats_samples", True)
        d["ycats_samples"] = self.pconfig.get("ycats_samples", True)
        return d

    def control_panel(self) -> str:
        """
        Heatmap-specific controls, only for the interactive version.
        """
        if self.flat:
            return ""

        # find min val across all datasets across all cols and rows
        minval = self.pconfig.get("min", None)
        if minval is None:
            minval = min(min(min(row) for row in dataset.data) for dataset in self.datasets)
        maxval = self.pconfig.get("max", None)
        if maxval is None:
            maxval = max(max(max(row) for row in dataset.data) for dataset in self.datasets)
        return f"""
        <div class="btn-group hc_switch_group">
            <button type="button" class="mqc_heatmap_sortHighlight btn btn-default btn-sm" data-target="#{self.id}" disabled="disabled">
                <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
            </button>
        </div>
        <div class="mqc_hcplot_range_sliders">
            <div>
                <label for="{self.id}_range_slider_min_txt">Min:</label>
                <input id="{self.id}_range_slider_min_txt" type="number" class="form-control" 
                    value="{minval}" data-target="{self.id}" data-minmax="min" min="{minval}" max="{maxval}" />
                <input id="{self.id}_range_slider_min" type="range" 
                    value="{minval}" data-target="{self.id}" data-minmax="min" min="{minval}" max="{maxval}" step="any" />
            </div>
            <div>
                <label for="{self.id}_range_slider_max_txt">Max:</label>
                <input id="{self.id}_range_slider_max_txt" type="number" class="form-control" 
                    value="{maxval}" data-target="{self.id}" data-minmax="max" min="{minval}" max="{maxval}" />
                <input id="{self.id}_range_slider_max" type="range" 
                    value="{maxval}" data-target="{self.id}" data-minmax="max" min="{minval}" max="{maxval}" step="any" />
            </div>
        </div>
        """

    def create_figure(
        self,
        layout: go.Layout,
        dataset: Dataset,
        is_log=False,
        is_pct=False,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        data: List[List[ElemT]] = dataset.data
        return go.Figure(
            data=go.Heatmap(
                z=data,
                x=self.xcats,
                y=self.ycats,
            ),
            layout=layout,
        )

    def save_data_file(self, dataset: Dataset) -> None:
        pass
