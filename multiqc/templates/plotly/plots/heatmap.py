import dataclasses
import logging
from typing import Dict, List, Union, Optional
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
from multiqc.utils import util_functions

logger = logging.getLogger(__name__)


ElemT = Union[str, float, int]


def plot(
    rows: List[List[ElemT]],
    xcats: List[str],
    ycats: List[str],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param rows: One dataset. A dataset is a list of rows of values
    :param xcats: Labels for X axis
    :param ycats: Labels for Y axis
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = HeatmapPlot(pconfig, rows, xcats, ycats)

    from multiqc.utils import report

    return p.add_to_report(report)


class HeatmapPlot(Plot):
    @dataclasses.dataclass
    class Dataset(BaseDataset):
        rows: List[List[ElemT]]
        xcats: List[str]
        ycats: List[str]
        square: bool = True

        @staticmethod
        def create(
            dataset: BaseDataset,
            pconfig: Dict,
            rows: List[List[ElemT]],
            xcats: List[str],
            ycats: List[str],
        ) -> "HeatmapPlot.Dataset":
            dataset = HeatmapPlot.Dataset(
                **dataset.__dict__,
                rows=rows,
                xcats=xcats,
                ycats=ycats,
                square=pconfig.get("square", True),  # Keep heatmap cells square
            )

            # Determining the size of the plot to reasonably display data without cluttering it too much.
            # For flat plots, we try to make the image large enough to display all samples, but to a limit
            # For interactive plots, we set a lower default height, as it will possible to resize the plot
            num_rows = len(ycats)
            num_cols = len(xcats)
            MIN_SIZE = 300
            MAX_HEIGHT = 900 if dataset.plot.flat else 500  # smaller number for interactive, as it's resizable
            MAX_WIDTH = 900  # default interactive width can be bigger

            font_size = 12
            # Making sure all ticks are displayed
            x_nticks = num_cols
            y_nticks = num_rows

            # Number of samples to the desired size in pixel one sample will take on a screen
            def n_elements_to_size(n: int):
                if n >= 50:
                    return 11
                if n >= 30:
                    return 18
                if n >= 20:
                    return 26
                if n >= 15:
                    return 30

            x_px_per_elem = n_elements_to_size(num_cols) or MIN_SIZE / num_cols
            y_px_per_elem = n_elements_to_size(num_rows) or MIN_SIZE / num_rows
            px_per_elem = min(x_px_per_elem, y_px_per_elem)
            if px_per_elem and px_per_elem <= 18:
                font_size = 10
            if px_per_elem and px_per_elem <= 13:
                font_size = 8

            height = num_rows * px_per_elem
            width = num_cols * px_per_elem

            if width < MIN_SIZE and height < MIN_SIZE:
                logger.debug(f"Resizing from {width}x{height} to fit the minimum size {MIN_SIZE}")
                px_per_elem = MIN_SIZE / max(num_rows, num_cols)
                width = num_cols * px_per_elem
                height = num_rows * px_per_elem
                logger.debug(f"Resized to {width}x{height}")
            if height > MAX_HEIGHT or width > MAX_WIDTH:
                logger.debug(f"Resizing from {width}x{height} to fit the maximum size {MAX_WIDTH}x{MAX_HEIGHT}")
                y_nticks = None  # allow skipping ticks to avoid making the font even smaller
                x_nticks = None  # allow skipping ticks to avoid making the font even smaller
                px_per_elem = min(MAX_WIDTH / num_cols, MAX_HEIGHT / num_rows)
                width = num_cols * px_per_elem
                height = num_rows * px_per_elem
                logger.debug(f"Resized to {width}x{height}")

            logger.debug(
                f"Heatmap size: {width}x{height}, px per element: {px_per_elem}, font: {font_size}px, xticks: {x_nticks}, yticks: {y_nticks}"
            )
            dataset.layout.xaxis.tickangle = 45
            dataset.layout.font.size = font_size
            dataset.layout.xaxis.nticks = x_nticks
            dataset.layout.yaxis.nticks = y_nticks
            dataset.layout.height = pconfig.get("height") or (200 + height)
            dataset.layout.width = pconfig.get("width") or (200 + width) if dataset.square else None
            dataset.layout.xaxis.automargin = True  # to make sure there is enough space for ticks labels
            dataset.layout.yaxis.automargin = True  # to make sure there is enough space for ticks labels
            dataset.layout.xaxis.showgrid = False
            dataset.layout.yaxis.showgrid = False
            dataset.layout.yaxis.autorange = "reversed"  # to make sure the first sample is at the top
            dataset.layout.yaxis.ticklabelposition = "outside right"

            colstops = pconfig.get("colstops", None)
            if colstops:
                # A list of 2-element lists where the first element is the
                # normalized color level value (starting at 0 and ending at 1),
                # and the second item is a valid color string.
                try:
                    colorscale = [[float(x), color] for x, color in colstops]
                except ValueError:
                    colorscale = None
                else:
                    # normalise the stop to a range from 0 to 1
                    minval = min(x for x, _ in colorscale)
                    maxval = max(x for x, _ in colorscale)
                    rng = maxval - minval
                    colorscale = [[(x - minval) / rng, color] for x, color in colorscale]
            else:
                # default colstops
                colorscale = [
                    [0, "#313695"],
                    [0.1, "#4575b4"],
                    [0.2, "#74add1"],
                    [0.3, "#abd9e9"],
                    [0.4, "#e0f3f8"],
                    [0.5, "#ffffbf"],
                    [0.6, "#fee090"],
                    [0.7, "#fdae61"],
                    [0.8, "#f46d43"],
                    [0.9, "#d73027"],
                    [1, "#a50026"],
                ]
            dataset.trace_params = {
                "colorscale": colorscale,
                "reversescale": pconfig.get("reverseColors", False),
                "showscale": pconfig.get("legend", True),
                "zmin": pconfig.get("min", None),
                "zmax": pconfig.get("max", None),
            }

            # Enable datalabels if there are less than 20x20 cells, unless heatmap_config.datalabels is set explicitly
            if pconfig.get("datalabels") is None and num_rows * num_cols < 400:
                dataset.trace_params["texttemplate"] = "%{z:.2f}"

            return dataset

        def create_figure(
            self,
            layout: Optional[go.Layout] = None,
            is_log=False,
            is_pct=False,
        ) -> go.Figure:
            """
            Create a Plotly figure for a dataset
            """
            return go.Figure(
                data=go.Heatmap(
                    z=self.rows,
                    x=self.xcats,
                    y=self.ycats,
                    **self.trace_params,
                ),
                layout=layout or self.layout,
            )

    def __init__(
        self,
        pconfig: Dict,
        rows: List[List[ElemT]],
        xcats: List[str],
        ycats: List[str],
    ):
        super().__init__(PlotType.HEATMAP, pconfig, n_datasets=1)

        # Extend each dataset object with a list of samples
        self.datasets: List[HeatmapPlot.Dataset] = [
            HeatmapPlot.Dataset.create(
                self.datasets[0],
                pconfig=pconfig,
                rows=rows,
                xcats=xcats,
                ycats=ycats,
            )
        ]

        self.min = self.pconfig.get("min", None)
        if self.min is None:
            for dataset in self.datasets:
                for row in dataset.rows:
                    for val in row:
                        if val is not None:
                            assert isinstance(val, (int, float))
                            self.min = val if self.min is None else min(self.min, val)
        self.max = self.pconfig.get("max", None)
        if self.max is None:
            for dataset in self.datasets:
                for row in dataset.rows:
                    for val in row:
                        if val is not None:
                            self.max = val if self.max is None else max(self.max, val)

    def dump_for_javascript(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().dump_for_javascript()
        d["xcats_samples"] = self.pconfig.get("xcats_samples", True)
        d["ycats_samples"] = self.pconfig.get("ycats_samples", True)
        return d

    def buttons(self) -> List[str]:
        """
        Heatmap-specific controls, only for the interactive version.
        """
        buttons = super().buttons()

        if not self.flat:
            # find min val across all datasets across all cols and rows
            buttons.append(
                f"""
            <div class="btn-group hc_switch_group">
                <button type="button" class="mqc_heatmap_sortHighlight btn btn-default btn-sm" data-target="#{self.id}" disabled="disabled">
                    <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
                </button>
            </div>
            <div class="mqc_hcplot_range_sliders">
                <div>
                    <label for="{self.id}_range_slider_min_txt">Min:</label>
                    <input id="{self.id}_range_slider_min_txt" type="number" class="form-control" 
                        value="{self.min}" data-target="{self.id}" data-minmax="min" min="{self.min}" max="{self.max}" />
                    <input id="{self.id}_range_slider_min" type="range" 
                        value="{self.min}" data-target="{self.id}" data-minmax="min" min="{self.min}" max="{self.max}" step="any" />
                </div>
                <div>
                    <label for="{self.id}_range_slider_max_txt">Max:</label>
                    <input id="{self.id}_range_slider_max_txt" type="number" class="form-control" 
                        value="{self.max}" data-target="{self.id}" data-minmax="max" min="{self.min}" max="{self.max}" />
                    <input id="{self.id}_range_slider_max" type="range" 
                        value="{self.max}" data-target="{self.id}" data-minmax="max" min="{self.min}" max="{self.max}" step="any" />
                </div>
            </div>
            """
            )
        return buttons

    def save_data_file(self, dataset: Dataset) -> None:
        data = [
            ["."] + dataset.xcats,
        ]
        for ycat, row in zip(dataset.ycats, dataset.rows):
            data.append([ycat] + row)

        util_functions.write_data_file(data, dataset.uid)