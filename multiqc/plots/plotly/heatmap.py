import dataclasses
import logging
from typing import Dict, List, Union, Optional
import plotly.graph_objects as go

from multiqc.plots.plotly.plot import Plot, PlotType, BaseDataset, split_long_string
from multiqc.utils import util_functions

logger = logging.getLogger(__name__)


ElemT = Union[str, float, int]


def plot(
    rows: Union[List[List[ElemT]], Dict[str, Dict[str, ElemT]]],
    pconfig: Dict,
    xcats: Optional[List[str]] = None,
    ycats: Optional[List[str]] = None,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param rows: One dataset. A dataset is a list of rows of values
    :param xcats: Labels for X axis
    :param ycats: Labels for Y axis
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = HeatmapPlot(rows, pconfig, xcats, ycats)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class Dataset(BaseDataset):
    rows: List[List[ElemT]]
    xcats: List[str]
    ycats: List[str]

    @staticmethod
    def create(
        dataset: BaseDataset,
        rows: Union[List[List[ElemT]], Dict[str, Dict[str, ElemT]]],
        xcats: Optional[List[str]] = None,
        ycats: Optional[List[str]] = None,
    ) -> "Dataset":
        if isinstance(rows, dict):
            # Convert dict to a list of lists
            if not ycats:
                ycats = list(rows.keys())
            if not xcats:
                xcats = []
                for y, value_by_x in rows.items():
                    for x, value in value_by_x.items():
                        if x not in xcats:
                            xcats.append(x)
            rows = [[rows.get(y, {}).get(x) for x in xcats] for y in ycats]

        dataset = Dataset(
            **dataset.__dict__,
            rows=rows,
            xcats=xcats,
            ycats=ycats,
        )
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


class HeatmapPlot(Plot):
    def __init__(
        self,
        rows: Union[List[List[ElemT]], Dict[str, Dict[str, ElemT]]],
        pconfig: Dict,
        xcats: Optional[List[str]],
        ycats: Optional[List[str]],
    ):
        super().__init__(PlotType.HEATMAP, pconfig, n_datasets=1)

        if isinstance(rows, list):
            if ycats and not isinstance(ycats, list):
                raise ValueError(
                    f"Heatmap plot {self.id}: ycats must be passed as a list when the input data is a 2d list. "
                    f"The order of that list should match the order of the rows in the input data."
                )
            if xcats and not isinstance(xcats, list):
                raise ValueError(
                    f"Heatmap plot {self.id}: xcats must be passed as a list when the input data is a 2d list. "
                    f"The order of that list should match the order of the columns in the input data."
                )

        self.layout.update(
            yaxis=dict(
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            showlegend=pconfig.get("legend", True),
        )

        self.square = pconfig.get("square", True)  # Keep heatmap cells square

        # Extend each dataset object with a list of samples
        self.datasets: List[Dataset] = [
            Dataset.create(
                self.datasets[0],
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
                        if val is not None and isinstance(val, (int, float)):
                            self.min = val if self.min is None else min(self.min, val)
        self.max = self.pconfig.get("max", None)
        if self.max is None:
            for dataset in self.datasets:
                for row in dataset.rows:
                    for val in row:
                        if val is not None and isinstance(val, (int, float)):
                            self.max = val if self.max is None else max(self.max, val)

        # Determining the size of the plot to reasonably display data without cluttering it too much.
        # For flat plots, we try to make the image large enough to display all samples, but to a limit
        # For interactive plots, we set a lower default height, as it will possible to resize the plot
        num_rows = len(self.datasets[0].ycats)
        num_cols = len(self.datasets[0].xcats)
        MAX_HEIGHT = 900 if self.flat else 500  # smaller number for interactive, as it's resizable
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
                return 60
            if n >= 3:
                return 80
            if n >= 2:
                return 100
            return 120

        x_px_per_elem = n_elements_to_size(num_cols)
        y_px_per_elem = n_elements_to_size(num_rows)
        min_px_per_elem = min(x_px_per_elem, y_px_per_elem)
        if self.square:
            x_px_per_elem = y_px_per_elem = min_px_per_elem

        width = pconfig.get("width") or int(num_cols * x_px_per_elem)
        height = pconfig.get("height") or int(num_rows * y_px_per_elem)

        if not self.square and width < MAX_WIDTH and x_px_per_elem < 40:  # can fit more columns on the screen
            logger.debug(f"Resizing width from {width} to {MAX_WIDTH} to fit horizontal column text on the screen")
            width = MAX_WIDTH
            x_px_per_elem = width / num_cols

        if height > MAX_HEIGHT or width > MAX_WIDTH:
            logger.debug(f"Resizing from {width}x{height} to fit the maximum size {MAX_WIDTH}x{MAX_HEIGHT}")
            if self.square:
                px_per_elem = min(MAX_WIDTH / num_cols, MAX_HEIGHT / num_rows)
                width = height = int(num_rows * px_per_elem)
            else:
                x_px_per_elem = MAX_WIDTH / num_cols
                y_px_per_elem = MAX_HEIGHT / num_rows
                width = int(num_cols * x_px_per_elem)
                height = int(num_rows * y_px_per_elem)

        logger.debug(f"Heatmap size: {width}x{height}, px per element: {x_px_per_elem:.2f}x{y_px_per_elem:.2f}")

        # For not very large datasets, making sure all ticks are displayed:
        if y_px_per_elem > 12:
            self.layout.yaxis.tickmode = "array"
            self.layout.yaxis.tickvals = list(range(num_rows))
            self.layout.yaxis.ticktext = ycats
        if x_px_per_elem > 18:
            self.layout.xaxis.tickmode = "array"
            self.layout.xaxis.tickvals = list(range(num_cols))
            self.layout.xaxis.ticktext = xcats
        if pconfig.get("angled_xticks", True) is False and x_px_per_elem >= 40:
            # Break up the horizontal ticks by whitespace to make them fit better vertically:
            self.layout.xaxis.ticktext = ["<br>".join(split_long_string(cat, 10)) for cat in xcats]
            # And leave x ticks horizontal:
            self.layout.xaxis.tickangle = 0
        else:
            # Rotate x-ticks to fit more of them on screen
            self.layout.xaxis.tickangle = 45

        self.layout.height = 200 + height
        self.layout.width = (250 + width) if self.square else None

        self.layout.xaxis.showgrid = False
        self.layout.yaxis.showgrid = False
        self.layout.yaxis.autorange = "reversed"  # to make sure the first sample is at the top
        self.layout.yaxis.ticklabelposition = "outside right"

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

        decimal_places = pconfig.get("tt_decimals", 2)

        xlab = pconfig.get("xlab", "x")
        ylab = pconfig.get("ylab", "y")
        zlab = pconfig.get("zlab", "z")
        hovertemplate = f"{xlab}: %{{x}}<br>{ylab}: %{{y}}<br>{zlab}: %{{z}}<extra></extra>"

        for ds in self.datasets:
            ds.trace_params = {
                "colorscale": colorscale,
                "reversescale": pconfig.get("reverseColors", False),
                "showscale": pconfig.get("legend", True),
                "zmin": pconfig.get("min", None),
                "zmax": pconfig.get("max", None),
                "hovertemplate": hovertemplate,
            }
            # Enable datalabels if there are less than 20x20 cells, unless heatmap_config.datalabels is set explicitly
            if pconfig.get("datalabels") is None and num_rows * num_cols < 400:
                ds.trace_params["texttemplate"] = "%{z:." + str(decimal_places) + "f}"

    def dump_for_javascript(self) -> Dict:
        """Serialise the plot data to pick up in JavaScript"""
        d = super().dump_for_javascript()
        d["xcats_samples"] = self.pconfig.get("xcats_samples", True)
        d["ycats_samples"] = self.pconfig.get("ycats_samples", True)
        d["square"] = self.square
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

    def save_data_file(self, dataset: Dataset) -> None:
        data = [
            ["."] + dataset.xcats,
        ]
        for ycat, row in zip(dataset.ycats, dataset.rows):
            data.append([ycat] + row)

        util_functions.write_data_file(data, dataset.uid)
