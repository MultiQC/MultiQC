"""MultiQC functions to plot a scatter plot"""

import copy
import logging
from collections import defaultdict
import stat
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple, Union, cast

import numpy as np
import pandas as pd
from plotly import graph_objects as go  # type: ignore

from multiqc import report
from multiqc.core.plot_data_store import parse_value
from multiqc.plots.plot import BaseDataset, NormalizedPlotInputData, PConfig, Plot, PlotType, plot_anchor
from multiqc.types import Anchor, SampleName

logger = logging.getLogger(__name__)


class ScatterConfig(PConfig):
    categories: Optional[List[str]] = None
    extra_series: Union[Dict[str, Any], List[Dict[str, Any]], List[List[Dict[str, Any]]], None] = None
    marker_size: Optional[int] = None
    marker_line_width: Optional[int] = None
    color: Optional[str] = None
    opacity: Optional[float] = None

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

    def to_df(self) -> pd.DataFrame:
        """
        Convert the scatter plot data to a pandas DataFrame for storage and reloading.
        """
        records = []

        # Serialize datasets, samplenames, and points
        for ds_idx, dataset in enumerate(self.datasets):
            dataset_label = self.extract_dataset_label(ds_idx)

            for sample_name, points in dataset.items():
                if not isinstance(points, list):
                    points = [points]

                for point_idx, point in enumerate(points):
                    record = {
                        "anchor": self.anchor,
                        "dataset_idx": ds_idx,
                        "dataset_label": dataset_label,
                        "sample_name": sample_name,
                        "point_idx": point_idx,
                    }

                    # Add all point data
                    for key, val in point.items():
                        if key == "x":
                            record["x"] = str(val)
                            record["x_type"] = type(val).__name__
                        if key == "y":
                            record["y"] = str(val)
                            record["y_type"] = type(val).__name__
                        else:
                            record[f"point_{key}"] = str(val)
                    records.append(record)

        df = pd.DataFrame(records)

        # Add config data as additional columns
        config_dict = self.pconfig.model_dump()
        for key, value in config_dict.items():
            # Only serialize primitive types directly
            if isinstance(value, (str, int, float, bool)) or value is None:
                df[f"config_{key}"] = value

        # Add anchor information
        df["anchor"] = str(self.anchor)

        return df

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
        )

    @classmethod
    def merge(
        cls, old_data: "ScatterNormalizedInputData", new_data: "ScatterNormalizedInputData"
    ) -> "ScatterNormalizedInputData":
        """
        Merge normalized data from old run and new run, matching by data labels when available
        """
        # Create dictionaries to map data_labels to datasets
        old_datasets_by_label = {}
        new_datasets_by_label = {}

        # Create a list for merged datasets
        merged_datasets = []

        # If data_labels exist, use them as keys for matching datasets
        if old_data.pconfig.data_labels and new_data.pconfig.data_labels:
            # First build mappings from label IDs to datasets
            for i, dl in enumerate(old_data.pconfig.data_labels):
                if i < len(old_data.datasets):
                    label_id = dl.get("name", f"dataset-{i}") if isinstance(dl, dict) else dl
                    old_datasets_by_label[label_id] = old_data.datasets[i]

            for i, dl in enumerate(new_data.pconfig.data_labels):
                if i < len(new_data.datasets):
                    label_id = dl.get("name", f"dataset-{i}") if isinstance(dl, dict) else dl
                    new_datasets_by_label[label_id] = new_data.datasets[i]

            # Create a dictionary to store merged data_labels
            merged_data_labels = []
            data_labels_by_id = {}

            # Store all data labels by their ID
            for i, dl in enumerate(old_data.pconfig.data_labels):
                if i < len(old_data.datasets):
                    label_id = dl.get("name", f"dataset-{i}") if isinstance(dl, dict) else dl
                    data_labels_by_id[label_id] = dl

            for i, dl in enumerate(new_data.pconfig.data_labels):
                if i < len(new_data.datasets):
                    label_id = dl.get("name", f"dataset-{i}") if isinstance(dl, dict) else dl
                    # New data labels override old ones with the same ID
                    data_labels_by_id[label_id] = dl

            # First process datasets that exist in both old and new data
            for label_id, old_ds in old_datasets_by_label.items():
                if label_id in new_datasets_by_label:
                    new_ds = new_datasets_by_label[label_id]
                    # Merge datasets - for scatter plots, combine all points
                    # Create a sample-to-points mapping from old dataset
                    sample_to_points = defaultdict(list)
                    for sample_name, points in old_ds.items():
                        if isinstance(points, list):
                            sample_to_points[sample_name].extend(points)
                        else:
                            sample_to_points[sample_name].append(points)

                    # Add points from new dataset
                    for sample_name, points in new_ds.items():
                        if isinstance(points, list):
                            sample_to_points[sample_name].extend(points)
                        else:
                            sample_to_points[sample_name].append(points)

                    # Convert back to a regular dict
                    merged_ds = dict(sample_to_points)
                    merged_datasets.append(merged_ds)
                    merged_data_labels.append(data_labels_by_id[label_id])

                    # Mark as processed
                    new_datasets_by_label.pop(label_id)
                else:
                    # Only in old data
                    merged_datasets.append(old_ds)
                    merged_data_labels.append(data_labels_by_id[label_id])

            # Then add datasets that only exist in new data
            for label_id, new_ds in new_datasets_by_label.items():
                merged_datasets.append(new_ds)
                merged_data_labels.append(data_labels_by_id[label_id])

            # Create a new pconfig with merged data_labels
            merged_pconf = ScatterConfig(**new_data.pconfig.model_dump())
            merged_pconf.data_labels = merged_data_labels

            return ScatterNormalizedInputData(
                plot_type=PlotType.SCATTER,
                anchor=new_data.anchor,
                datasets=merged_datasets,
                pconfig=merged_pconf,
            )
        else:
            # If no data labels, preserve old behavior (return new data)
            # But combine extra_series if present
            if old_data.pconfig.extra_series and new_data.pconfig.extra_series is None:
                merged_pconf = ScatterConfig(**new_data.pconfig.model_dump())
                merged_pconf.extra_series = old_data.pconfig.extra_series
                return ScatterNormalizedInputData(
                    plot_type=PlotType.SCATTER,
                    anchor=new_data.anchor,
                    datasets=new_data.datasets,
                    pconfig=merged_pconf,
                )
            return new_data

    @classmethod
    def from_df(
        cls, df: pd.DataFrame, pconfig: Union[Dict, ScatterConfig], anchor: Anchor
    ) -> "ScatterNormalizedInputData":
        """
        Create a ScatterNormalizedInputData object from a pandas DataFrame.
        """
        # Filter out rows related to this anchor if there are multiple
        if "anchor" in df.columns:
            df = df[df["anchor"] == str(anchor)]

        pconf = (
            pconfig
            if isinstance(pconfig, ScatterConfig)
            else cast(ScatterConfig, ScatterConfig.from_pconfig_dict(pconfig))
        )

        if df.empty:
            # Return empty data if no valid rows found
            return ScatterNormalizedInputData(
                plot_type=PlotType.SCATTER,
                anchor=anchor,
                datasets=[],
                pconfig=pconf,
            )

        # Extract config information that might have been serialized
        config_cols = [col for col in df.columns if col.startswith("config_")]
        config_data = {}
        for col in config_cols:
            key = col[7:]  # Remove "config_" prefix
            if df[col].nunique() == 1:
                config_data[key] = df[col].iloc[0]

        pconf = cast(ScatterConfig, ScatterConfig.from_pconfig_dict({**config_data, **pconf.model_dump()}))

        # Reconstruct datasets
        dataset_indices = sorted(df["dataset_idx"].unique())
        datasets = []

        for dataset_idx in dataset_indices:
            dataset_df = df[df["dataset_idx"] == dataset_idx]
            dataset = {}

            # Group by sample name
            for sample_name in dataset_df["sample_name"].unique():
                sample_df = dataset_df[dataset_df["sample_name"] == sample_name]

                # Extract points
                points = []
                for _, row in sample_df.iterrows():
                    point = {}
                    for col in row.index:
                        point["x"] = parse_value(row["x"], row["x_type"])
                        point["y"] = parse_value(row["y"], row["y_type"])
                        if col.startswith("point_"):
                            key = col[6:]  # Remove "point_" prefix
                            point[key] = row[col]
                    points.append(point)

                dataset[sample_name] = points

            datasets.append(dataset)

        cls.dataset_labels_from_df(df, pconf)

        return cls(
            anchor=anchor,
            plot_type=PlotType.SCATTER,
            datasets=datasets,
            pconfig=pconf,
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
            color = el.get("color")
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
            n_samples_per_dataset=[len(x) for x in points_lists],
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

                    for k in ["color", "opacity", "marker_size", "marker_line_width"]:
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
        inputs.save()
        return plot
