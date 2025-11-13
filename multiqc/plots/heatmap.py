"""MultiQC functions to plot a heatmap"""

import logging
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union, cast

import numpy as np
import plotly.graph_objects as go  # type: ignore
import polars as pl
from pydantic import Field

from multiqc import report
from multiqc.core.plot_data_store import parse_value
from multiqc.plots.plot import (
    BaseDataset,
    NormalizedPlotInputData,
    PConfig,
    Plot,
    PlotType,
    plot_anchor,
    split_long_string,
)
from multiqc.types import Anchor, SampleName
from multiqc.utils.util_functions import scipy_hierarchy_leaves_list, scipy_hierarchy_linkage, scipy_pdist

logger = logging.getLogger(__name__)


# Define element types for the heatmap
ElemT = Union[str, float, int, None]


class HeatmapConfig(PConfig):
    """Configuration for a heatmap plot"""

    xlab: str = "x"
    ylab: str = "y"
    zlab: str = "z"
    min: Union[float, int, None] = None
    max: Union[float, int, None] = None
    xcats_samples: bool = True
    ycats_samples: bool = True
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
    cluster_rows: bool = True
    cluster_cols: bool = True
    cluster_method: str = "complete"  # linkage method: single, complete, average, weighted, etc.
    cluster_switch_clustered_active: bool = False

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("heatmap",), **data)


class HeatmapNormalizedInputData(NormalizedPlotInputData):
    """
    Represents normalized input data for a heatmap plot.
    """

    rows: List[List[ElemT]]
    xcats: List[Union[str, int]]
    ycats: List[Union[str, int]]
    pconfig: HeatmapConfig

    def is_empty(self) -> bool:
        return len(self.rows) == 0 or all(len(row) == 0 for row in self.rows)

    def to_df(self) -> pl.DataFrame:
        """
        Convert the heatmap data to a polars DataFrame for storage and reloading.
        """
        # Create a DataFrame with row indices, column indices, and values
        records = []
        for i, row in enumerate(self.rows):
            y_cat = str(self.ycats[i]) if i < len(self.ycats) else str(i)
            for j, val in enumerate(row):
                x_cat = str(self.xcats[j]) if j < len(self.xcats) else str(j)
                records.append(
                    {
                        "row_idx": i,
                        "col_idx": j,
                        "row_cat": y_cat,
                        "col_cat": x_cat,
                        # values can be be different types (int, float, str...), especially across
                        # plots. parquet requires values of the same type. so we cast them to str
                        "z_val": str(val),
                        "z_val_type": type(val).__name__,
                    }
                )

        df = pl.DataFrame(records)
        return self.finalize_df(df)

    @classmethod
    def from_df(
        cls, df: pl.DataFrame, pconfig: Union[Dict, HeatmapConfig], anchor: Anchor
    ) -> "HeatmapNormalizedInputData":
        """
        Create a HeatmapNormalizedInputData object from a polars DataFrame.
        """
        # Filter out rows related to this anchor if there are multiple
        if "anchor" in df.columns:
            df = df.filter(pl.col("anchor") == str(anchor))

        if df.is_empty():
            # Return empty data if no valid rows found
            pconf = (
                pconfig
                if isinstance(pconfig, HeatmapConfig)
                else cast(HeatmapConfig, HeatmapConfig.from_pconfig_dict(pconfig))
            )
            return cls(
                anchor=anchor,
                rows=[],
                xcats=[],
                ycats=[],
                pconfig=pconf,
                plot_type=PlotType.HEATMAP,
                creation_date=cls.creation_date_from_df(df),
            )

        pconf = cast(HeatmapConfig, HeatmapConfig.from_df(df))

        # Extract and rebuild the 2D matrix
        max_row_idx = int(df.select(pl.col("row_idx").max()).item())
        max_col_idx = int(df.select(pl.col("col_idx").max()).item())

        # Reconstruct xcats and ycats
        xcats = []
        ycats = []

        # Get unique row and column categories in the order of indices
        row_cats_df = df.select(["row_idx", "row_cat"]).unique().sort("row_idx")
        col_cats_df = df.select(["col_idx", "col_cat"]).unique().sort("col_idx")

        ycats = row_cats_df.select("row_cat").to_series().to_list()
        xcats = col_cats_df.select("col_cat").to_series().to_list()

        # Create empty matrix with appropriate type annotations
        rows: List[List[ElemT]] = []
        for _ in range(max_row_idx + 1):
            row: List[ElemT] = [None] * (max_col_idx + 1)  # type: ignore
            rows.append(row)

        # Fill matrix with values
        for r in df.iter_rows(named=True):
            row_idx = int(r["row_idx"])
            col_idx = int(r["col_idx"])
            val = parse_value(r["z_val"], r["z_val_type"])
            rows[row_idx][col_idx] = val

        return cls(
            anchor=anchor,
            rows=rows,
            xcats=xcats,
            ycats=ycats,
            pconfig=pconf,
            plot_type=PlotType.HEATMAP,
            creation_date=cls.creation_date_from_df(df),
        )

    @classmethod
    def merge(
        cls, old_data: "HeatmapNormalizedInputData", new_data: "HeatmapNormalizedInputData"
    ) -> "HeatmapNormalizedInputData":
        """
        Merge two HeatmapNormalizedInputData objects.
        This method is called when merging data from multiple runs.

        The merge takes the union of all categories - if old_data has A and B,
        and new_data has A and C, the result will have A, B, and C.

        For any overlapping x-y pairs, the value from new_data (latest) is used.
        """
        # Check for empty data cases
        if old_data.is_empty():
            return new_data
        if new_data.is_empty():
            return old_data

        # Create a union of x and y categories
        all_xcats = list(dict.fromkeys(old_data.xcats + new_data.xcats))
        all_ycats = list(dict.fromkeys(old_data.ycats + new_data.ycats))

        # Create mappings from categories to indices
        xcat_to_idx = {str(cat): i for i, cat in enumerate(all_xcats)}
        ycat_to_idx = {str(cat): i for i, cat in enumerate(all_ycats)}

        # Initialize a matrix with None values
        merged_rows: List[List[ElemT]] = []
        for _ in range(len(all_ycats)):
            row: List[ElemT] = [None] * len(all_xcats)  # type: ignore
            merged_rows.append(row)

        # Helper function to fill the matrix from a data source
        def fill_matrix(data, is_newer=False):
            for y_idx, row in enumerate(data.rows):
                if y_idx >= len(data.ycats):
                    continue  # Skip extra rows
                y_cat = str(data.ycats[y_idx])
                if y_cat not in ycat_to_idx:
                    continue  # Skip unknown categories

                merged_y_idx = ycat_to_idx[y_cat]

                for x_idx, val in enumerate(row):
                    if x_idx >= len(data.xcats):
                        continue  # Skip extra columns
                    x_cat = str(data.xcats[x_idx])
                    if x_cat not in xcat_to_idx:
                        continue  # Skip unknown categories

                    merged_x_idx = xcat_to_idx[x_cat]

                    # For newer data or missing values, always update
                    if is_newer or merged_rows[merged_y_idx][merged_x_idx] is None:
                        merged_rows[merged_y_idx][merged_x_idx] = val

        # Fill matrix with values - old data first, then new data to override duplicates
        fill_matrix(old_data, is_newer=False)
        fill_matrix(new_data, is_newer=True)

        # Use the newer config (from new_data)
        return HeatmapNormalizedInputData(
            anchor=new_data.anchor,
            rows=merged_rows,
            xcats=all_xcats,
            ycats=all_ycats,
            pconfig=new_data.pconfig,
            plot_type=PlotType.HEATMAP,
            creation_date=report.creation_date,
        )

    @staticmethod
    def create(
        data: Union[Sequence[Sequence[ElemT]], Mapping[Union[str, int], Mapping[Union[str, int], ElemT]]],
        xcats: Optional[Sequence[Union[str, int]]] = None,
        ycats: Optional[Sequence[Union[str, int]]] = None,
        pconfig: Union[Dict[str, Any], HeatmapConfig, None] = None,
    ) -> "HeatmapNormalizedInputData":
        pconf = cast(HeatmapConfig, HeatmapConfig.from_pconfig_dict(pconfig))

        rows: List[List[ElemT]]
        if isinstance(data, dict):
            # Re-key the dict to be strings
            rows_str: Dict[str, Dict[str, ElemT]] = {
                str(y): {str(x): value for x, value in value_by_x.items()} for y, value_by_x in data.items()
            }

            # Convert dict to a list of lists
            if not ycats:
                ycats = list(rows_str.keys())
            if not xcats:
                xcats = []
                for _, value_by_x in rows_str.items():
                    for x, _ in value_by_x.items():
                        if x not in xcats:
                            xcats.append(x)
            rows = [[rows_str.get(str(y), {}).get(str(x)) for x in xcats] for y in ycats]
        else:
            rows = cast(List[List[ElemT]], data)
            if ycats is None:
                ycats = xcats
            if xcats is None:
                xcats = list(range(len(rows[0])))
            if ycats is None:
                ycats = list(range(len(rows)))

            if ycats and not isinstance(ycats, list):
                raise ValueError(
                    "Heatmap plot: ycats must be passed as a list when the input data is a 2d list. "
                    "The order of that list should match the order of the rows in the input data."
                )
            if xcats and not isinstance(xcats, list):
                raise ValueError(
                    "Heatmap plot: xcats must be passed as a list when the input data is a 2d list. "
                    "The order of that list should match the order of the columns in the input data."
                )

        return HeatmapNormalizedInputData(
            anchor=plot_anchor(pconf),
            plot_type=PlotType.HEATMAP,
            rows=rows,
            xcats=list(xcats),
            ycats=list(ycats),
            pconfig=pconf,
            creation_date=report.creation_date,
        )


def plot(
    data: Union[Sequence[Sequence[ElemT]], Mapping[Union[str, int], Mapping[Union[str, int], ElemT]]],
    xcats: Optional[Sequence[Union[str, int]]] = None,
    ycats: Optional[Sequence[Union[str, int]]] = None,
    pconfig: Union[Dict[str, Any], HeatmapConfig, None] = None,
) -> Union["HeatmapPlot", str, None]:
    """
    Plot a 2D heatmap.
    :param data: List of lists, each a representing a row of values; or a dict of dicts
    :param xcats: Labels for x-axis
    :param ycats: Labels for y-axis. Defaults to same as x.
    :param pconfig: optional dict with config key:value pairs.
    :return: HTML and JS, ready to be inserted into the page
    """
    inputs: HeatmapNormalizedInputData = HeatmapNormalizedInputData.create(data, xcats, ycats, pconfig)
    inputs = HeatmapNormalizedInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    return HeatmapPlot.from_inputs(inputs)


def _cluster_data(
    rows: List[List[ElemT]], cluster_rows: bool = True, cluster_cols: bool = True, method: str = "complete"
) -> Tuple[List[List[ElemT]], List[int], List[int]]:
    """Cluster the heatmap data and return clustered data with new indices"""
    row_idx = list(range(len(rows)))
    col_idx = list(range(len(rows[0])))

    data_array = np.array([[0.0 if x is None else float(x) for x in row] for row in rows])

    if cluster_rows and len(rows) > 1:
        try:
            row_dist = scipy_pdist(data_array)
            row_linkage = scipy_hierarchy_linkage(row_dist, method=method)
            row_idx = scipy_hierarchy_leaves_list(row_linkage)
            data_array = data_array[row_idx]
        except Exception as e:
            logger.warning(f"Row clustering failed: {str(e)}")

    if cluster_cols and len(rows[0]) > 1:
        try:
            col_dist = scipy_pdist(data_array.T)
            col_linkage = scipy_hierarchy_linkage(col_dist, method=method)
            col_idx = scipy_hierarchy_leaves_list(col_linkage)
            data_array = data_array[:, col_idx]
        except Exception as e:
            logger.warning(f"Column clustering failed: {str(e)}")

    return cast(List[List[ElemT]], data_array.tolist()), row_idx, col_idx


class Dataset(BaseDataset):
    rows: List[List[ElemT]]
    rows_clustered: Optional[List[List[ElemT]]] = None
    xcats: Sequence[str]
    ycats: Sequence[str]
    xcats_clustered: Optional[Sequence[str]] = None
    ycats_clustered: Optional[Sequence[str]] = None
    xcats_samples: bool = True
    ycats_samples: bool = True

    def sample_names(self) -> List[SampleName]:
        snames: List[SampleName] = []
        if self.xcats_samples:
            snames.extend(SampleName(cat) for cat in self.xcats)
        if self.ycats_samples:
            snames.extend(SampleName(cat) for cat in self.ycats)
        return snames

    @staticmethod
    def create(
        dataset: BaseDataset,
        rows: List[List[ElemT]],
        xcats: Sequence[Union[str, int]],
        ycats: Sequence[Union[str, int]],
        cluster_rows: bool = True,
        cluster_cols: bool = True,
        cluster_method: str = "complete",
        xcats_samples: bool = True,
        ycats_samples: bool = True,
    ) -> "Dataset":
        rows_clustered = None
        xcats_clustered = None
        ycats_clustered = None

        if cluster_rows or cluster_cols:
            try:
                clustered_rows, row_idx, col_idx = _cluster_data(rows, cluster_rows, cluster_cols, cluster_method)
                rows_clustered = clustered_rows
                xcats_clustered = [xcats[i] for i in col_idx] if cluster_cols else xcats
                ycats_clustered = [ycats[i] for i in row_idx] if cluster_rows else ycats
            except Exception as e:
                logger.warning(f"Clustering failed: {str(e)}")

        dataset = Dataset(
            **dataset.__dict__,
            rows=rows,
            rows_clustered=rows_clustered,
            xcats=[str(x) for x in xcats],
            ycats=[str(y) for y in ycats],
            xcats_clustered=[str(x) for x in xcats_clustered] if xcats_clustered else None,
            ycats_clustered=[str(y) for y in ycats_clustered] if ycats_clustered else None,
            xcats_samples=xcats_samples,
            ycats_samples=ycats_samples,
        )
        return dataset

    def create_figure(
        self,
        layout: Optional[go.Layout] = None,
        is_log: bool = False,
        is_pct: bool = False,
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
        if self.xcats:
            row += self.xcats
        data: List[List[ElemT]] = []
        data.append(row)
        for i, row in enumerate(self.rows):
            new_row: List[ElemT] = []
            if self.ycats:
                ycat = self.ycats[i]
                new_row.append(ycat)
            new_row += row
            data.append(new_row)

        report.write_data_file(data, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: HeatmapConfig, keep_hidden: bool = True) -> str:
        """
        Format heatmap as a markdown table
        """
        prompt = ""
        if self.xcats:
            if self.ycats:
                prompt = "|"
                if pconfig.ycats_samples:
                    prompt += "Sample"
            xcats = self.xcats
            if pconfig.xcats_samples:
                xcats = [report.anonymize_sample_name(cat) for cat in xcats]
            prompt += "|" + "|".join(xcats) + "|\n"
            if self.ycats:
                prompt += "|---"
            prompt += "|" + "|".join("---" for _ in self.xcats) + "|\n"
        for i, row in enumerate(self.rows):
            if self.ycats:
                ycat = self.ycats[i]
                if pconfig.ycats_samples:
                    ycat = report.anonymize_sample_name(ycat)
                    prompt += "|" + ycat
            prompt += "|" + "|".join(str(x) for x in row) + "|\n"
        return prompt


class HeatmapPlot(Plot[Dataset, HeatmapConfig]):
    datasets: List[Dataset]
    xcats_samples: bool
    ycats_samples: bool
    min: Optional[float] = None
    max: Optional[float] = None
    cluster_switch_clustered_active: bool = False

    def sample_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        if self.xcats_samples:
            for ds in self.datasets:
                if ds.xcats:
                    names.extend(SampleName(cat) for cat in ds.xcats)
        if self.ycats_samples:
            for ds in self.datasets:
                if ds.ycats:
                    names.extend(SampleName(cat) for cat in ds.ycats)
        return names

    @staticmethod
    def from_inputs(inputs: HeatmapNormalizedInputData) -> Union["HeatmapPlot", str, None]:
        plot = HeatmapPlot.create(
            rows=inputs.rows,
            pconfig=inputs.pconfig,
            anchor=inputs.anchor,
            xcats=list(inputs.xcats),
            ycats=list(inputs.ycats),
        )
        inputs.save_to_parquet()
        return plot

    @staticmethod
    def create(
        rows: List[List[ElemT]],
        pconfig: HeatmapConfig,
        anchor: Anchor,
        xcats: List[Union[str, int]],
        ycats: List[Union[str, int]],
    ) -> "HeatmapPlot":
        max_n_rows = 0
        max_n_cols = 0
        if rows:
            max_n_rows = len(rows)
            if len(rows[0]) > 0:
                max_n_cols = max(max_n_cols, len(rows[0]))

        n_samples = max_n_cols if pconfig.xcats_samples else (max_n_rows if pconfig.ycats_samples else 0)

        model: Plot[Dataset, HeatmapConfig] = Plot.initialize(
            plot_type=PlotType.HEATMAP,
            pconfig=pconfig,
            anchor=anchor,
            n_series_per_dataset=[max_n_rows],
            n_samples_per_dataset=[n_samples],
            defer_render_if_large=False,  # We hide samples on large heatmaps, so no need to defer render
            flat_if_very_large=True,  # However, the data is still embedded into the HTML, and we don't want the report size to inflate
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
        xcats_samples: bool = pconfig.xcats_samples
        ycats_samples: bool = pconfig.ycats_samples

        # Extend each dataset object with a list of samples
        model.datasets = [
            Dataset.create(
                model.datasets[0],
                rows=rows,
                xcats=xcats,
                ycats=ycats,
                cluster_rows=pconfig.cluster_rows,
                cluster_cols=pconfig.cluster_cols,
                cluster_method=pconfig.cluster_method,
                xcats_samples=xcats_samples,
                ycats_samples=ycats_samples,
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
        num_rows = len(model.datasets[0].ycats) if model.datasets[0].ycats else len(model.datasets[0].rows)
        num_cols = len(model.datasets[0].xcats) if model.datasets[0].xcats else len(model.datasets[0].rows[0])
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
            cluster_switch_clustered_active=pconfig.cluster_switch_clustered_active,
        )

    def buttons(self, flat: bool, module_anchor: Anchor, section_anchor: Anchor) -> List[str]:
        """
        Heatmap-specific controls, only for the interactive version.

        """
        buttons = super().buttons(flat=flat, module_anchor=module_anchor, section_anchor=section_anchor)
        if any(ds.rows_clustered for ds in self.datasets):
            # find min val across all datasets across all cols and rows
            buttons.append(
                f"""
                <div class="btn-group" role="group">
                    <button
                        type="button"
                        class="btn btn-default btn-sm {"" if self.pconfig.cluster_switch_clustered_active else "active"}"
                        data-action="unclustered"
                        data-plot-anchor="{self.anchor}"
                    >
                        Sorted by sample
                    </button>
                    <button
                        type="button"
                        class="btn btn-default btn-sm {"active" if self.pconfig.cluster_switch_clustered_active else ""}"
                        data-action="clustered"
                        data-plot-anchor="{self.anchor}"
                    >
                        Clustered
                    </button>
                </div>
                """
            )

        return buttons

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.xlab:
            result += f"X axis: {self.pconfig.xlab}\n"
        if self.pconfig.ylab:
            result += f"Y axis: {self.pconfig.ylab}\n"
        if self.pconfig.zlab:
            result += f"Z axis: {self.pconfig.zlab}\n"
        return result
