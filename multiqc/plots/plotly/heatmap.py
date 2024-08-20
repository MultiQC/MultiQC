import logging
from typing import Dict, List, Union, Optional, Tuple
import plotly.graph_objects as go  # type: ignore
from pydantic import Field

from multiqc.plots.plotly.plot import (
    PlotType,
    BaseDataset,
    split_long_string,
    Plot,
    PConfig,
)
from multiqc import report

logger = logging.getLogger(__name__)


ElemT = Union[str, float, int, None]


class HeatmapConfig(PConfig):
    """Configuration for a heatmap plot"""

    xlab: str = "x"
    ylab: str = "y"
    zlab: str = "z"
    min: Union[float, int, None] = None
    max: Union[float, int, None] = None
    xcats_samples: bool = False
    ycats_samples: bool = False
    square: bool = True
    colstops: List[List] = []
    reverseColors: bool = Field(False, deprecated="reverse_colors")
    reverse_colors: bool = False
    decimalPlaces: int = Field(2, deprecated="tt_decimals")
    tt_decimals: int = 2
    legend: bool = True
    datalabels: Optional[bool] = Field(None, deprecated="display_values")
    display_values: Optional[bool] = None
    angled_xticks: bool = True


def plot(
    rows: Union[List[List[ElemT]], Dict[Union[str, int], Dict[Union[str, int], ElemT]]],
    pconfig: HeatmapConfig,
    xcats: Optional[List[Union[str, int]]] = None,
    ycats: Optional[List[Union[str, int]]] = None,
) -> "HeatmapPlot":
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param rows: One dataset. A dataset is a list of rows of values
    :param xcats: Labels for X axis
    :param ycats: Labels for Y axis
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    return HeatmapPlot.create(rows, pconfig, xcats, ycats)


class Dataset(BaseDataset):
    rows: List[List[ElemT]]
    xcats: List[str]
    ycats: List[str]

    @staticmethod
    def create(
        dataset: BaseDataset,
        rows: Union[List[List[ElemT]], Dict[Union[str, int], Dict[Union[str, int], ElemT]]],
        xcats: Optional[List[Union[str, int]]] = None,
        ycats: Optional[List[Union[str, int]]] = None,
    ) -> "Dataset":
        if isinstance(rows, dict):
            # Re-key the dict to be strings
            rows_str: Dict[str, Dict[str, ElemT]] = {
                str(y): {str(x): value for x, value in value_by_x.items()} for y, value_by_x in rows.items()
            }

            # Convert dict to a list of lists
            if not ycats:
                ycats = list(rows_str.keys())
            if not xcats:
                xcats = []
                for y, value_by_x in rows_str.items():
                    for x, value in value_by_x.items():
                        if x not in xcats:
                            xcats.append(x)
            rows = [[rows_str.get(str(y), {}).get(str(x)) for x in xcats] for y in ycats]

        dataset = Dataset(
            **dataset.__dict__,
            rows=rows,
            xcats=[str(x) for x in xcats] if xcats else None,
            ycats=[str(y) for y in ycats] if ycats else None,
        )
        return dataset

    def create_figure(
        self,
        layout: Optional[go.Layout] = None,
        is_log=False,
        is_pct=False,
        **kwargs,
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

    def save_data_file(self) -> None:
        row: List[ElemT] = ["."]
        row += self.xcats
        data: List[List[ElemT]] = []
        data.append(row)
        for ycat, row in zip(self.ycats, self.rows):
            new_row: List[ElemT] = [ycat]
            new_row += row
            data.append(new_row)

        report.write_data_file(data, self.uid)


class HeatmapPlot(Plot):
    datasets: List[Dataset]
    xcats_samples: bool
    ycats_samples: bool
    min: Optional[float] = None
    max: Optional[float] = None

    @staticmethod
    def create(
        rows: Union[List[List[ElemT]], Dict[Union[str, int], Dict[Union[str, int], ElemT]]],
        pconfig: HeatmapConfig,
        xcats: Optional[List[Union[str, int]]],
        ycats: Optional[List[Union[str, int]]],
    ) -> "HeatmapPlot":
        max_n_samples = 0
        if rows:
            max_n_samples = len(rows)
            if isinstance(rows, list):
                max_n_samples = max(max_n_samples, len(rows[0]))
            else:
                max_n_samples = max(max_n_samples, len(rows[next(iter(rows))]))

        model = Plot.initialize(
            plot_type=PlotType.HEATMAP,
            pconfig=pconfig,
            n_samples_per_dataset=[max_n_samples],
            defer_render_if_large=False,  # We hide samples on large heatmaps, so no need to defer render
            flat_if_very_large=True,  # However, the data is still embedded into the HTML, and we don't want the report size to inflate
        )

        if isinstance(rows, list):
            if ycats and not isinstance(ycats, list):
                raise ValueError(
                    f"Heatmap plot {model.id}: ycats must be passed as a list when the input data is a 2d list. "
                    f"The order of that list should match the order of the rows in the input data."
                )
            if xcats and not isinstance(xcats, list):
                raise ValueError(
                    f"Heatmap plot {model.id}: xcats must be passed as a list when the input data is a 2d list. "
                    f"The order of that list should match the order of the columns in the input data."
                )

        model.layout.update(
            yaxis=dict(
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            showlegend=pconfig.legend,
        )

        model.square = pconfig.square  # Keep heatmap cells square
        xcats_samples = pconfig.xcats_samples
        ycats_samples = pconfig.ycats_samples

        # Extend each dataset object with a list of samples
        model.datasets = [
            Dataset.create(
                model.datasets[0],
                rows=rows,
                xcats=xcats,
                ycats=ycats,
            )
        ]

        minval = model.pconfig.min
        if minval is None:
            for dataset in model.datasets:
                for row in dataset.rows:
                    for val in row:
                        if val is not None and isinstance(val, (int, float)):
                            minval = val if minval is None else min(minval, val)  # type: ignore
        maxval = model.pconfig.max
        if maxval is None:
            for dataset in model.datasets:
                for row in dataset.rows:
                    for val in row:
                        if val is not None and isinstance(val, (int, float)):
                            maxval = val if maxval is None else max(maxval, val)  # type: ignore

        # Determining the size of the plot to reasonably display data without cluttering it too much.
        # For flat plots, we try to make the image large enough to display all samples, but to a limit
        # For interactive plots, we set a lower default height, as it will possible to resize the plot
        num_rows = len(model.datasets[0].ycats)
        num_cols = len(model.datasets[0].xcats)
        MAX_HEIGHT = 900 if model.flat else 500  # smaller number for interactive, as it's resizable
        MAX_WIDTH = 900  # default interactive width can be bigger

        # Number of samples to the desired size in pixel one sample will take on a screen
        def n_elements_to_size(n: int):
            if n >= 50:
                return 11
            if n >= 30:
                return 18
            if n >= 20:
                return 26
            if n >= 15:
                return 33
            if n >= 10:
                return 45
            if n >= 5:
                return 52
            if n >= 3:
                return 60
            if n >= 2:
                return 80
            return 100

        x_px_per_elem = n_elements_to_size(num_cols)
        y_px_per_elem = n_elements_to_size(num_rows)
        min_px_per_elem = min(x_px_per_elem, y_px_per_elem)
        if model.square:
            x_px_per_elem = y_px_per_elem = min_px_per_elem

        width = pconfig.width or int(num_cols * x_px_per_elem)
        height = pconfig.height or int(num_rows * y_px_per_elem)

        if not model.square and width < MAX_WIDTH and x_px_per_elem < 40:  # can fit more columns on the screen
            # logger.debug(f"Resizing width from {width} to {MAX_WIDTH} to fit horizontal column text on the screen")
            width = MAX_WIDTH
            x_px_per_elem = width / num_cols

        if height > MAX_HEIGHT or width > MAX_WIDTH:
            # logger.debug(f"Resizing from {width}x{height} to fit the maximum size {MAX_WIDTH}x{MAX_HEIGHT}")
            if model.square:
                px_per_elem = min(MAX_WIDTH / num_cols, MAX_HEIGHT / num_rows)
                width = height = int(num_rows * px_per_elem)
            else:
                x_px_per_elem = MAX_WIDTH / num_cols
                y_px_per_elem = MAX_HEIGHT / num_rows
                width = int(num_cols * x_px_per_elem)
                height = int(num_rows * y_px_per_elem)

        # logger.debug(f"Heatmap size: {width}x{height}, px per element: {x_px_per_elem:.2f}x{y_px_per_elem:.2f}")

        # For not very large datasets, making sure all ticks are displayed:
        if y_px_per_elem > 12:
            model.layout.yaxis.tickmode = "array"
            model.layout.yaxis.tickvals = list(range(num_rows))
            model.layout.yaxis.ticktext = ycats
        if x_px_per_elem > 18:
            model.layout.xaxis.tickmode = "array"
            model.layout.xaxis.tickvals = list(range(num_cols))
            model.layout.xaxis.ticktext = xcats
        if not pconfig.angled_xticks and x_px_per_elem >= 40 and xcats:
            # Break up the horizontal ticks by whitespace to make them fit better vertically:
            model.layout.xaxis.ticktext = ["<br>".join(split_long_string(str(cat), 10)) for cat in xcats]
            # And leave x ticks horizontal:
            model.layout.xaxis.tickangle = 0
        else:
            # Rotate x-ticks to fit more of them on screen
            model.layout.xaxis.tickangle = 45

        model.layout.height = 200 + height
        model.layout.width = (250 + width) if model.square else None

        model.layout.xaxis.showgrid = False
        model.layout.yaxis.showgrid = False
        model.layout.yaxis.autorange = "reversed"  # to make sure the first sample is at the top
        model.layout.yaxis.ticklabelposition = "outside right"

        colorscale: List[Tuple[float, str]] = []
        if pconfig.colstops:
            # A list of 2-element lists where the first element is the
            # normalized color level value (starting at 0 and ending at 1),
            # and the second item is a valid color string.
            try:
                colorscale = [(float(x), color) for x, color in pconfig.colstops]
            except ValueError:
                pass
            else:
                # normalise the stop to a range from 0 to 1
                minval = min(x for x, _ in colorscale)
                maxval = max(x for x, _ in colorscale)
                rng = maxval - minval
                colorscale = [((x - minval) / rng, color) for x, color in colorscale]
        else:
            # default colstops
            colorscale = [
                (0, "#313695"),
                (0.1, "#4575b4"),
                (0.2, "#74add1"),
                (0.3, "#abd9e9"),
                (0.4, "#e0f3f8"),
                (0.5, "#ffffbf"),
                (0.6, "#fee090"),
                (0.7, "#fdae61"),
                (0.8, "#f46d43"),
                (0.9, "#d73027"),
                (1, "#a50026"),
            ]

        xlab = pconfig.xlab
        ylab = pconfig.ylab
        zlab = pconfig.zlab
        hovertemplate = f"{xlab}: %{{x}}<br>{ylab}: %{{y}}<br>{zlab}: %{{z}}<extra></extra>"

        for ds in model.datasets:
            ds.trace_params = {
                "colorscale": colorscale,
                "reversescale": pconfig.reverse_colors,
                "showscale": pconfig.legend,
                "zmin": pconfig.min,
                "zmax": pconfig.max,
                "hovertemplate": hovertemplate,
            }
            # Enable labels if there are less than 20x20 cells, unless display_values is set explicitly
            if pconfig.display_values is True or pconfig.display_values is None and num_rows * num_cols < 400:
                ds.trace_params["texttemplate"] = "%{z:." + str(pconfig.tt_decimals) + "f}"

            # We only want to use heatmap labels for tooltips, but not on the axis
            ds.layout["xaxis"]["title"] = None
            ds.layout["yaxis"]["title"] = None

        return HeatmapPlot(
            **model.__dict__,
            xcats_samples=xcats_samples,
            ycats_samples=ycats_samples,
            min=minval,
            max=maxval,
        )

    def buttons(self, flat: bool) -> List[str]:
        """
        Heatmap-specific controls, only for the interactive version.
        """
        buttons = super().buttons(flat=flat)

        if not flat:
            # find min val across all datasets across all cols and rows
            buttons.append(
                f"""
            <div class="mqc_hcplot_range_sliders">
                <div>
                    <label for="{self.id}_range_slider_min_txt">Min:</label>
                    <input id="{self.id}_range_slider_min_txt" type="number" class="form-control" 
                        value="{self.min}" data-target="{self.id}" data-minmax="min" min="{self.min}" max="{self.max}" />
                    <input id="{self.id}_range_slider_min" type="range" 
                        value="{self.min}" data-target="{self.id}" data-minmax="min" min="{self.min}" max="{self.max}" step="any" />
                </div>
                <div style="margin-left: 30px;">
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
