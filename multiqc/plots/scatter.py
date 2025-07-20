"""MultiQC functions to plot a scatter plot"""

import copy
import json
import logging
import math
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple, Union, cast

import numpy as np
import polars as pl
from plotly import graph_objects as go  # type: ignore

from multiqc import report
from multiqc.core.plot_data_store import parse_value
from multiqc.plots.plot import BaseDataset, NormalizedPlotInputData, PConfig, Plot, PlotType, plot_anchor
from multiqc.types import Anchor, SampleName
from multiqc.utils import mqc_colour

logger = logging.getLogger(__name__)


class ScatterConfig(PConfig):
    categories: Optional[List[str]] = None
    extra_series: Union[Dict[str, Any], List[Dict[str, Any]], List[List[Dict[str, Any]]], None] = None
    marker_size: Optional[int] = None
    marker_line_width: Optional[int] = None
    color: Optional[str] = None
    opacity: Optional[float] = None
    marker_symbol: Optional[str] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("scatterplot",), **data)


# {'color': 'rgb(211,211,211,0.05)', 'name': 'background: EUR', 'x': -0.294, 'y': -1.527}
ValueT = Union[str, float, int]
PointT = Dict[str, ValueT]


class ScatterNormalizedInputData(NormalizedPlotInputData):
    datasets: List[Dict[str, Any]]
    pconfig: ScatterConfig

    def is_empty(self) -> bool:
        return len(self.datasets) == 0 or all(len(ds) == 0 for ds in self.datasets)

    def to_df(self) -> pl.DataFrame:
        """
        Convert the scatter plot data to a polars DataFrame for storage and reloading.
        """
        records = []

        # Serialize datasets, samplenames, and points
        for ds_idx, dataset in enumerate(self.datasets):
            for sample_name, points in dataset.items():
                if not isinstance(points, list):
                    points = [points]

                for point_idx, point in enumerate(points):
                    record = {
                        "dataset_idx": ds_idx,
                        "data_label": json.dumps(self.pconfig.data_labels[ds_idx]) if self.pconfig.data_labels else "",
                        "sample": sample_name,
                        "point_idx": point_idx,
                    }

                    # Add all point data
                    for key, val in point.items():
                        if key == "x":
                            # Convert NaN values to string marker for safe serialization
                            x_val = "__NAN__MARKER__" if isinstance(val, float) and math.isnan(val) else str(val)
                            record["x"] = x_val
                            record["x_type"] = type(val).__name__
                        elif key == "y":
                            # Convert NaN values to string marker for safe serialization
                            y_val = "__NAN__MARKER__" if isinstance(val, float) and math.isnan(val) else str(val)
                            record["y"] = y_val
                            record["y_type"] = type(val).__name__
                        else:
                            # For other values, also handle NaN
                            point_val = "__NAN__MARKER__" if isinstance(val, float) and math.isnan(val) else str(val)
                            record[f"point_{key}"] = point_val
                    records.append(record)

        df = pl.DataFrame(
            records,
            schema_overrides={
                "data_label": pl.Utf8,
                "sample": pl.Utf8,
            },
        )
        return self.finalize_df(df)

    @staticmethod
    def create(
        data: Union[Dict[str, Any], List[Dict[str, Any]]],
        pconfig: Union[Mapping[str, Any], ScatterConfig, None] = None,
    ) -> "ScatterNormalizedInputData":
        pconf: ScatterConfig = cast(ScatterConfig, ScatterConfig.from_pconfig_dict(pconfig))

        # Given one dataset - turn it into a list
        if not isinstance(data, list):
            data = [data]  # type: ignore

        return ScatterNormalizedInputData(
            anchor=plot_anchor(pconf),
            plot_type=PlotType.SCATTER,
            datasets=data,
            pconfig=pconf,
            creation_date=report.creation_date,
        )

    @classmethod
    def merge(
        cls, old_data: "ScatterNormalizedInputData", new_data: "ScatterNormalizedInputData"
    ) -> "ScatterNormalizedInputData":
        """
        Merge normalized data from old run and new run, matching by data labels when available
        """
        new_df = new_data.to_df()
        if new_df.is_empty():
            return old_data

        old_df = old_data.to_df()

        # If we have both old and new data, merge them
        merged_df = new_df
        if old_df is not None and not old_df.is_empty():
            # Get the list of samples that exist in both old and new data, for each dataset
            new_keys = new_df.select(["data_label", "sample"]).unique()

            # Keep only the rows in old_df whose (data_label, sample) pair
            # does *not* appear in new_keys. An anti-join is the cleanest way to express this.
            old_df_filtered = old_df.join(
                new_keys,
                on=["data_label", "sample"],
                how="anti",  # anti-join = “rows in left not matched in right”
            )

            # Combine the filtered old data with new data
            merged_df = pl.concat([old_df_filtered, new_df], how="diagonal")

        return cls.from_df(merged_df, new_data.pconfig, new_data.anchor)

    @classmethod
    def from_df(
        cls, df: pl.DataFrame, pconfig: Union[Dict, ScatterConfig], anchor: Anchor
    ) -> "ScatterNormalizedInputData":
        """
        Create a ScatterNormalizedInputData object from a polars DataFrame.
        """
        if df.is_empty():
            pconf = (
                pconfig
                if isinstance(pconfig, ScatterConfig)
                else cast(ScatterConfig, ScatterConfig.from_pconfig_dict(pconfig))
            )
            return cls(
                anchor=anchor,
                datasets=[],
                pconfig=pconf,
                plot_type=PlotType.SCATTER,
                creation_date=cls.creation_date_from_df(df),
            )

        pconf = cast(ScatterConfig, ScatterConfig.from_df(df))

        # Reconstruct datasets
        data_labels = []
        datasets = []
        dataset_indices = sorted(df.select("dataset_idx").unique().to_series().to_list())

        for dataset_idx in dataset_indices:
            dataset_df = df.filter(pl.col("dataset_idx") == dataset_idx)
            dataset = {}

            data_label = dataset_df.select("data_label").row(0)[0]
            data_labels.append(json.loads(data_label) if data_label else {})

            # Group by sample name
            sample_names = dataset_df.select("sample").unique().to_series().to_list()
            for sample_name in sample_names:
                sample_df = dataset_df.filter(pl.col("sample") == sample_name)

                # Extract points
                points = []
                for row in sample_df.iter_rows(named=True):
                    point = {}
                    point["x"] = parse_value(row["x"], row["x_type"])
                    point["y"] = parse_value(row["y"], row["y_type"])
                    # Add any additional point attributes
                    for col in row.keys():
                        if col.startswith("point_"):
                            key = col[6:]  # Remove "point_" prefix
                            point[key] = row[col]
                    points.append(point)

                dataset[sample_name] = points

            datasets.append(dataset)

        if any(d for d in data_labels if d):
            pconf.data_labels = data_labels
        return cls(
            anchor=anchor,
            plot_type=PlotType.SCATTER,
            datasets=datasets,
            pconfig=pconf,
            creation_date=cls.creation_date_from_df(df),
        )


def plot(
    data: Union[Dict[str, Any], List[Dict[str, Any]]],
    pconfig: Union[Mapping[str, Any], ScatterConfig, None],
) -> Union["ScatterPlot", str, None]:
    """
    Plot a scatter plot with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    inputs: ScatterNormalizedInputData = ScatterNormalizedInputData.create(data, pconfig)
    inputs = ScatterNormalizedInputData.merge_with_previous(inputs)

    if inputs.is_empty():
        return None

    return ScatterPlot.from_inputs(inputs)


class Dataset(BaseDataset):
    points: List[PointT]

    def sample_names(self) -> List[SampleName]:
        return [
            SampleName(point["name"]) for point in self.points if "name" in point and isinstance(point["name"], str)
        ]

    @staticmethod
    def create(
        dataset: BaseDataset,
        points: List[Dict[str, Any]],
        pconfig: ScatterConfig,
    ) -> "Dataset":
        dataset = Dataset(
            **dataset.__dict__,
            points=points,
        )

        dataset.trace_params.update(
            textfont=dict(size=8),
            marker=dict(
                size=10,
                line=dict(width=1),
                opacity=1,
                color="rgba(124, 181, 236, .5)",
                symbol="circle",
            ),
        )
        # if categories is provided, set them as x-axis ticks
        if pconfig.categories:
            dataset.layout["xaxis"]["tickmode"] = "array"
            dataset.layout["xaxis"]["tickvals"] = list(range(len(pconfig.categories)))
            dataset.layout["xaxis"]["ticktext"] = pconfig.categories
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)
        MAX_ANNOTATIONS = 10  # Maximum number of dots to be annotated directly on the plot
        n_annotated = len([el for el in self.points if "annotation" in el])
        if n_annotated < MAX_ANNOTATIONS:
            points = [(i, x) for i, x in enumerate(self.points) if x.get("annotate", True) is not False]
            x_values = np.array([x["x"] for (i, x) in points])
            y_values = np.array([x["y"] for (i, x) in points])

            x_are_numeric = np.issubdtype(x_values.dtype, np.number)
            y_are_numeric = np.issubdtype(y_values.dtype, np.number)

            if x_are_numeric and y_are_numeric:
                # Finding and marking outliers to only label them
                # 1. Calculate Z-scores
                x_std = np.std(x_values)
                y_std = np.std(y_values)
                if x_std == 0 and y_std == 0:
                    logger.warning(f"Scatter plot {self.plot_id}: all {len(points)} points have the same coordinates")
                    if len(points) == 1:  # Only single point - annotate it!
                        for _, point in points:
                            point["annotation"] = point["name"]
                else:
                    x_z_scores = np.abs((x_values - np.mean(x_values)) / x_std) if x_std else np.zeros_like(x_values)
                    y_z_scores = np.abs((y_values - np.mean(y_values)) / y_std) if y_std else np.zeros_like(y_values)
                    # 2. Find a reasonable threshold so there are not too many outliers
                    threshold = 1.0
                    while threshold <= 6.0:
                        n_outliers = np.count_nonzero((x_z_scores > threshold) | (y_z_scores > threshold))
                        # logger.debug(f"Scatter plot outlier threshold: {threshold:.2f}, outliers: {n_outliers}")
                        if n_annotated + n_outliers <= MAX_ANNOTATIONS:
                            break
                        # If there are too many outliers, we increase the threshold until we have less than 10
                        threshold += 0.2
                    # 3. Annotate outliers that pass the threshold
                    for (_, point), x_z_score, y_z_score in zip(points, x_z_scores, y_z_scores):
                        # Check if point is an outlier or if total points are less than 10
                        if x_z_score > threshold or y_z_score > threshold:
                            point["annotation"] = point["name"]
                n_annotated = len([point for point in self.points if "annotation" in point])

        # If there are few unique colors, we can additionally put a unique list into a legend
        # (even though some color might belong to many distinct names - we will just crop the list)
        names_by_legend_key: Dict[Tuple[Any, Any, Any, Any], Set[str]] = defaultdict(set)

        for el in self.points:
            legend_key = (el.get("color"), el.get("marker_size"), el.get("marker_line_width"), el.get("group"))
            name = el["name"]
            assert isinstance(name, str)
            names_by_legend_key[legend_key].add(name)
        layout.showlegend = True

        in_legend: Set[Tuple[Any, Any, Any, Any]] = set()
        for el in self.points:
            x = el["x"]
            name = el["name"]
            group = el.get("group")
            color = mqc_colour.color_to_rgb_string(cast(Optional[str], el.get("color")))
            annotation = el.get("annotation")

            show_in_legend = False
            if layout.showlegend and not el.get("hide_in_legend"):
                key = (color, el.get("marker_size"), el.get("marker_line_width"), group)
                if key not in in_legend:
                    in_legend.add(key)
                    names = sorted(names_by_legend_key[key])
                    label = ", ".join(names)
                    if group:
                        label = f"{group}: {label}"
                    MAX_WIDTH = 60
                    if len(label) > MAX_WIDTH:  # Crop long labels
                        label = label[:MAX_WIDTH] + "..."
                    show_in_legend = True
                    name = label

            params = copy.deepcopy(self.trace_params)
            marker = params.pop("marker")
            if color:
                marker["color"] = color

            if "marker_line_width" in el:
                marker["line"]["width"] = el["marker_line_width"]
            if "marker_size" in el:
                marker["size"] = el["marker_size"]
            if "marker_symbol" in el:
                marker["symbol"] = el["marker_symbol"]
            if "opacity" in el:
                marker["opacity"] = el["opacity"]

            if annotation:
                params["mode"] = "markers+text"

            if n_annotated > 0:  # Reduce opacity of the borders that clutter the annotations:
                marker["line"]["color"] = "rgba(0, 0, 0, .2)"

            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[el["y"]],
                    name=name,
                    text=[annotation or name],
                    showlegend=show_in_legend,
                    marker=marker,
                    **params,
                )
            )
        fig.layout.height += len(in_legend) * 5  # extra space for legend
        return fig

    def get_x_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        if not self.points:
            return None, None
        xmax, xmin = None, None
        for point in self.points:
            x = point["x"]
            if xmax is not None:
                xmax = max(xmax, x)
            else:
                xmax = x
            if xmin is not None:
                xmin = min(xmin, x)
            else:
                xmin = x
        return xmin, xmax

    def get_y_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        if not self.points:
            return None, None
        ymax, ymin = None, None
        for point in self.points:
            y = point["y"]
            if ymax is not None:
                ymax = max(ymax, y)
            else:
                ymax = y
            if ymin is not None:
                ymin = min(ymin, y)
            else:
                ymin = y
        return ymin, ymax

    def save_data_file(self) -> None:
        data = [
            {
                "Name": point["name"],
                "X": point["x"],
                "Y": point["y"],
            }
            for point in self.points
        ]
        report.write_data_file(data, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        xsuffix = self.layout.get("xaxis", {}).get("ticksuffix", "")
        ysuffix = self.layout.get("yaxis", {}).get("ticksuffix", "")

        prompt = ""
        prompt += "|Sample|X|Y|\n"
        prompt += "|---|---|---|\n"

        for point in self.points:
            pseudonym = report.anonymize_sample_name(cast(str, point["name"]))
            prompt += f"|{pseudonym}|{point['x']}{xsuffix}|{point['y']}{ysuffix}|\n"

        return prompt


class ScatterPlot(Plot[Dataset, ScatterConfig]):
    datasets: List[Dataset]

    @staticmethod
    def create(
        points_lists: List[List[PointT]],
        pconfig: ScatterConfig,
        anchor: Anchor,
    ) -> "ScatterPlot":
        model: Plot[Dataset, ScatterConfig] = Plot.initialize(
            plot_type=PlotType.SCATTER,
            pconfig=pconfig,
            anchor=anchor,
            n_series_per_dataset=[len(x) for x in points_lists],
            default_tt_label="<br><b>X</b>: %{x}<br><b>Y</b>: %{y}",
        )

        model.datasets = [Dataset.create(d, points, pconfig) for d, points in zip(model.datasets, points_lists)]

        model._set_x_bands_and_range(pconfig)
        model._set_y_bands_and_range(pconfig)

        # Make a tooltip always show on hover over nearest point on plot
        model.layout.hoverdistance = -1

        return ScatterPlot(**model.__dict__)

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.xlab:
            result += f"X axis: {self.pconfig.xlab}\n"
        if self.pconfig.ylab:
            result += f"Y axis: {self.pconfig.ylab}\n"
        if self.pconfig.categories:
            result += f"X categories: {', '.join(self.pconfig.categories)}\n"

        return result

    @staticmethod
    def from_inputs(inputs: ScatterNormalizedInputData) -> Union["ScatterPlot", str, None]:
        pconf = inputs.pconfig
        sample_names = []
        plotdata: List[List[Dict[str, Any]]] = list()
        for data_index, ds in enumerate(inputs.datasets):
            d: List[Dict[str, Any]] = list()
            for s_name in ds:
                sample_names.append(SampleName(s_name))
                # Ensure any overwriting conditionals from data_labels (e.g. ymax) are taken in consideration
                series_config: ScatterConfig = pconf.model_copy()
                if pconf.data_labels:
                    dl = pconf.data_labels[data_index]
                    if isinstance(dl, dict):
                        # if not a dict: only dataset name is provided
                        for k, v in dl.items():
                            if k in series_config.model_fields:
                                setattr(series_config, k, v)

                if not isinstance(ds[s_name], list):
                    ds[s_name] = [ds[s_name]]
                for point in ds[s_name]:
                    if point["x"] is not None:
                        if series_config.xmax is not None and float(point["x"]) > float(series_config.xmax):
                            continue
                        if series_config.xmin is not None and float(point["x"]) < float(series_config.xmin):
                            continue
                    if point["y"] is not None:
                        if series_config.ymax is not None and float(point["y"]) > float(series_config.ymax):
                            continue
                        if series_config.ymin is not None and float(point["y"]) < float(series_config.ymin):
                            continue
                    if "name" in point:
                        point["name"] = f"{s_name}: {point['name']}"
                    else:
                        point["name"] = s_name

                    for k in ["color", "opacity", "marker_size", "marker_line_width", "marker_symbol"]:
                        if k not in point:
                            v = getattr(series_config, k)
                            if v is not None:
                                if isinstance(v, dict) and s_name in v:
                                    point[k] = v[s_name]
                                else:
                                    point[k] = v
                    d.append(point)

            plotdata.append(d)

        if pconf.square:
            if pconf.ymax is None and pconf.xmax is None:
                # Find the max value
                max_val = 0.0
                for d in plotdata:
                    for s in d:
                        max_val = max(max_val, s["x"], s["y"])
                max_val = 1.02 * float(max_val)  # add 2% padding
                pconf.xmax = pconf.xmax if pconf.xmax is not None else max_val
                pconf.ymax = pconf.ymax if pconf.ymax is not None else max_val

        # Add extra annotation data series
        # noinspection PyBroadException
        try:
            if pconf.extra_series:
                extra_series: List[List[Dict[str, Any]]] = []
                if isinstance(pconf.extra_series, dict):
                    extra_series = [[pconf.extra_series]]
                elif isinstance(pconf.extra_series[0], dict):
                    extra_series = [cast(List[Dict[str, Any]], [pconf.extra_series])]
                else:
                    extra_series = cast(List[List[Dict[str, Any]]], pconf.extra_series)
                for i, es in enumerate(extra_series):
                    for s in es:
                        plotdata[i].append(s)

        except Exception:
            pass

        plot = ScatterPlot.create(
            points_lists=plotdata,
            pconfig=pconf,
            anchor=inputs.anchor,
        )
        inputs.save_to_parquet()
        return plot
