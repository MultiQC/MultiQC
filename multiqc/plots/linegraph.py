"""MultiQC functions to plot a linegraph"""

import io
import json
import logging
import math
import os
import random
from typing import Any, Dict, Generic, List, Literal, Mapping, Optional, Sequence, Tuple, Type, TypeVar, Union, cast

import plotly.graph_objects as go  # type: ignore
import polars as pl
from natsort import natsorted
from pydantic import Field

from multiqc import config, report
from multiqc.core.plot_data_store import parse_value
from multiqc.plots.plot import (
    BaseDataset,
    NormalizedPlotInputData,
    PConfig,
    Plot,
    PlotType,
    convert_dash_style,
    plot_anchor,
)
from multiqc.types import Anchor, SampleName
from multiqc.utils import mqc_colour
from multiqc.utils.util_functions import update_dict
from multiqc.validation import ValidatedConfig, add_validation_warning

logger = logging.getLogger(__name__)


KeyT = TypeVar("KeyT", int, str, float)
ValT = TypeVar("ValT", int, str, float, None)
XToYDictT = Mapping[KeyT, ValT]
DatasetT = Mapping[Union[str, SampleName], XToYDictT[KeyT, ValT]]


class Marker(ValidatedConfig):
    symbol: Optional[str] = None
    color: Optional[str] = None
    line_color: Optional[str] = None
    fill_color: Optional[str] = None
    width: int = 1

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("Marker",), **data)


class Series(ValidatedConfig, Generic[KeyT, ValT]):
    name: str = Field(default_factory=lambda: f"series-{random.randint(1000000, 9999999)}")
    pairs: List[Tuple[KeyT, ValT]]
    color: Optional[str] = None
    width: int = 2
    dash: Optional[str] = None
    showlegend: bool = True
    marker: Optional[Marker] = None
    # Store additional trace parameters that should be passed to Plotly
    extra_trace_params: Dict[str, Any] = Field(default_factory=dict)

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        path_in_cfg = path_in_cfg or ("Series",)

        if "dashStyle" in data:
            add_validation_warning(path_in_cfg, "'dashStyle' field is deprecated. Please use 'dash' instead")
            data["dash"] = data.pop("dashStyle")

        tuples: List[Tuple[KeyT, ValT]] = []
        if "data" in data:
            add_validation_warning(path_in_cfg + ("data",), "'data' field is deprecated. Please use 'pairs' instead")
        for p in data.pop("data") if "data" in data else data.get("pairs", []):
            if isinstance(p, list):
                tuples.append(tuple(p))
            else:
                tuples.append(p)
        data["pairs"] = tuples

        # Extract extra trace parameters (fields not in the main model)
        main_fields = {"name", "pairs", "color", "width", "dash", "showlegend", "marker", "extra_trace_params"}
        extra_params = {k: v for k, v in data.items() if k not in main_fields}
        if extra_params:
            data["extra_trace_params"] = extra_params
            # Remove extra params from data to avoid validation errors
            for k in extra_params:
                data.pop(k)

        super().__init__(**data, path_in_cfg=path_in_cfg)

        if self.dash is not None:
            self.dash = convert_dash_style(self.dash, path_in_cfg=path_in_cfg + ("dash",))

    def get_x_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        xs = [x[0] for x in self.pairs]
        if len(xs) > 0:
            return min(xs), max(xs)  # type: ignore
        return None, None

    def get_y_range(self) -> Tuple[Optional[Any], Optional[Any]]:
        ys = [x[1] for x in self.pairs if x[1] is not None]
        if len(ys) > 0:
            return min(ys), max(ys)  # type: ignore
        return None, None


SeriesT = Union[Series, Dict[str, Any]]


AxisStr = Literal["xaxis", "yaxis"]


class LinePlotConfig(PConfig):
    xlab: Optional[str] = None
    ylab: Optional[str] = None
    categories: bool = False
    smooth_points: Optional[int] = 500
    smooth_points_sumcounts: Union[bool, List[bool], None] = None
    extra_series: Optional[Union[Series, List[Series], List[List[Series]]]] = None
    style: Optional[Literal["lines", "lines+markers"]] = None
    hide_empty: bool = Field(True)
    colors: Dict[str, str] = {}
    dash_styles: Dict[str, str] = {}
    hovertemplates: Dict[str, str] = {}
    legend_groups: Dict[str, str] = {}
    axis_controlled_by_switches: Optional[List[AxisStr]] = None

    @classmethod
    def parse_extra_series(
        cls,
        data: Union[SeriesT, List[SeriesT], List[List[SeriesT]]],
        path_in_cfg: Tuple[str, ...],
    ) -> Union[Series, List[Series], List[List[Series]]]:
        if isinstance(data, list):
            if isinstance(data[0], list):
                return [[Series(path_in_cfg=path_in_cfg, **d) if isinstance(d, dict) else d for d in ds] for ds in data]  # type: ignore
            return [Series(path_in_cfg=path_in_cfg, **d) if isinstance(d, dict) else d for d in data]  # type: ignore
        return Series(path_in_cfg=path_in_cfg, **data) if isinstance(data, dict) else data  # type: ignore

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("lineplot",), **data)


class Dataset(BaseDataset, Generic[KeyT, ValT]):
    lines: List[Series[KeyT, ValT]]

    def sample_names(self) -> List[SampleName]:
        return [SampleName(line.name) for line in self.lines]

    def get_x_range(self) -> Tuple[Optional[KeyT], Optional[KeyT]]:
        if not self.lines:
            return None, None
        xmax, xmin = None, None
        for line in self.lines:
            _xmin, _xmax = line.get_x_range()
            if _xmin is not None:
                xmin = min(xmin, _xmin) if xmin is not None else _xmin  # type: ignore
            if _xmax is not None:
                xmax = max(xmax, _xmax) if xmax is not None else _xmax  # type: ignore
        return xmin, xmax

    def get_y_range(self) -> Tuple[Optional[ValT], Optional[ValT]]:
        if not self.lines:
            return None, None
        ymax, ymin = None, None
        for line in self.lines:
            _ymin, _ymax = line.get_y_range()
            if _ymin is not None:
                ymin = min(ymin, _ymin) if ymin is not None else _ymin  # type: ignore
            if _ymax is not None:
                ymax = max(ymax, _ymax) if ymax is not None else _ymax  # type: ignore
        return ymin, ymax

    @staticmethod
    def create(
        base_dataset: BaseDataset,
        lines: List[Series[KeyT, ValT]],
        pconfig: LinePlotConfig,
    ) -> "Dataset[KeyT, ValT]":
        dataset: Dataset[KeyT, ValT] = Dataset(**base_dataset.model_dump(), lines=lines)

        # Prevent Plotly-JS from parsing strings as numbers
        if pconfig.categories or dataset.dconfig.get("categories"):
            dataset.layout["xaxis"]["type"] = "category"

        if pconfig.style is not None:
            mode = pconfig.style
        else:
            num_data_points = sum(len(x.pairs) for x in lines)
            if num_data_points < config.lineplot_number_of_points_to_hide_markers:
                mode = "lines+markers"
            else:
                mode = "lines"

        dataset.trace_params.update(
            mode=mode,
            line={"width": 2},
        )
        if mode == "lines+markers":
            dataset.trace_params.update(
                line={"width": 0.6},
                marker={"size": 5},
            )
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
        if layout.showlegend is True:
            # Extra space for legend
            if hasattr(layout, "height") and isinstance(layout.height, int):
                layout.height += len(self.lines) * 5

        fig = go.Figure(layout=layout)
        for series in self.lines:
            xs = [x[0] for x in series.pairs]
            ys = [x[1] for x in series.pairs]
            params: Dict[str, Any] = {
                "showlegend": series.showlegend,
                "line": {
                    "color": series.color,
                    "dash": series.dash,
                    "width": series.width,
                },
            }
            if series.marker:
                params["mode"] = "lines+markers"
                params["marker"] = {
                    "symbol": series.marker.symbol,
                    "color": series.marker.fill_color or series.marker.color or series.color,
                    "line": {
                        "width": series.marker.width,
                        "color": series.marker.line_color or series.marker.color or "black",
                    },
                }
            params = update_dict(params, self.trace_params, none_only=True)
            # Add extra trace parameters from series
            params = update_dict(params, series.extra_trace_params, none_only=True)
            if len(series.pairs) == 1:
                params["mode"] = "lines+markers"  # otherwise it's invisible

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=series.name,
                    text=[series.name] * len(xs),
                    **params,
                )
            )
        return fig

    def save_data_file(self) -> None:
        y_by_x_by_sample: Dict[str, Dict[Union[float, str], Any]] = dict()
        last_cats = None
        shared_cats = True
        for series in self.lines:
            y_by_x_by_sample[series.name] = dict()

            # Check to see if all categories are the same
            if len(series.pairs) > 0 and isinstance(series.pairs[0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in series.pairs]
                elif last_cats != [x[0] for x in series.pairs]:
                    shared_cats = False

            for i, x in enumerate(series.pairs):
                if isinstance(x, list):
                    y_by_x_by_sample[series.name][x[0]] = x[1]
                else:
                    try:
                        y_by_x_by_sample[series.name][self.dconfig["categories"][i]] = x
                    except (ValueError, KeyError, IndexError):
                        y_by_x_by_sample[series.name][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format in ["tsv", "csv"]:
            sep = "\t" if config.data_format == "tsv" else ","
            fout = ""
            for series in self.lines:
                fout += series.name + sep + "X" + sep + sep.join([str(x[0]) for x in series.pairs]) + "\n"
                fout += series.name + sep + "Y" + sep + sep.join([str(x[1]) for x in series.pairs]) + "\n"

            fn = f"{self.uid}.{config.data_format_extensions[config.data_format]}"
            fpath = os.path.join(report.data_tmp_dir(), fn)
            with io.open(fpath, "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            report.write_data_file(y_by_x_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        xsuffix = self.layout.get("xaxis", {}).get("ticksuffix", "")
        ysuffix = self.layout.get("yaxis", {}).get("ticksuffix", "")

        # Use pseudonyms for sample names if available
        pseudonyms = [report.anonymize_sample_name(series.name) for series in self.lines]

        # Create header with axis information and common suffixes
        result = "Samples: " + ", ".join(pseudonyms) + "\n\n"

        # If all y-values have the same suffix (like %), mention it in the header
        if ysuffix:
            result += f"Y values are in {ysuffix}\n\n"
        if xsuffix:
            result += f"X values are in {xsuffix}\n\n"

        for pseudonym, series in zip(pseudonyms, self.lines):
            # For other plots or fewer points, use the original format but without redundant suffixes
            pairs = [f"{self.fmt_value_for_llm(x[0])}: {self.fmt_value_for_llm(x[1])}" for x in series.pairs]

            result += f"{pseudonym} {', '.join(pairs)}\n\n"

        return result


class LinePlotNormalizedInputData(NormalizedPlotInputData[LinePlotConfig], Generic[KeyT, ValT]):
    """
    Represents normalized input data for a line plot.

    We want to be permissive with user input, e.g. allow one dataset or a list of datasets,
    allow optional categories, categories as strings or as dicts or as objects. We want
    to normalize the input data to dump it to parequet, before we use it to plot.
    """

    data: List[List[Series[KeyT, ValT]]]
    sample_names: List[SampleName]

    def is_empty(self) -> bool:
        return len(self.data) == 0 or all(len(ds) == 0 for ds in self.data)

    def to_df(self) -> pl.DataFrame:
        """
        Save plot data to a parquet file using a tabular representation that's
        optimized for cross-run analysis.

        Instead of serializing complex nested structures, we create a structured
        table with explicit columns for series data.
        """
        records = []
        # Create a record for each data point in each series
        for ds_idx, dataset in enumerate(self.data):
            data_label = json.dumps(self.pconfig.data_labels[ds_idx]) if self.pconfig.data_labels else ""
            for series in dataset:
                # Extract series properties once per series, not per data point
                series_props = {k: v for k, v in series.model_dump().items() if k not in ["pairs", "name"]}
                sample_name = series.name
                for x, y in series.pairs:
                    # Convert NaN values to string marker for safe serialization
                    x_val = "__NAN__MARKER__" if isinstance(x, float) and math.isnan(x) else str(x)
                    y_val = "__NAN__MARKER__" if isinstance(y, float) and math.isnan(y) else str(y)

                    record = {
                        "dataset_idx": ds_idx,
                        "data_label": data_label,
                        "sample": sample_name,
                        # values can be be different types (int, float, str...), especially across
                        # plots. parquet requires values of the same type. so we cast them to str
                        "x_val": x_val,
                        "y_val": y_val,
                        "x_val_type": type(x).__name__,
                        "y_val_type": type(y).__name__,
                        "series": series_props,
                    }
                    records.append(record)

        # Create DataFrame from records
        df = pl.DataFrame(
            records,
            schema_overrides={
                "data_label": pl.Utf8,
                "sample": pl.Utf8,
                "x_val": pl.Utf8,
                "y_val": pl.Utf8,
                "x_val_type": pl.Utf8,
                "y_val_type": pl.Utf8,
            },
        )
        return self.finalize_df(df)

    @classmethod
    def from_df(
        cls, df: pl.DataFrame, pconfig: Union[Dict, LinePlotConfig], anchor: Anchor
    ) -> "LinePlotNormalizedInputData[KeyT, ValT]":
        pconf: LinePlotConfig
        if cls.df_is_empty(df):
            pconf = (
                pconfig
                if isinstance(pconfig, LinePlotConfig)
                else cast(LinePlotConfig, LinePlotConfig.from_pconfig_dict(pconfig))
            )
            return cls(
                anchor=anchor,
                plot_type=PlotType.LINE,
                data=[],
                pconfig=pconf,
                sample_names=[],
                creation_date=cls.creation_date_from_df(df),
            )
        pconf = cast(LinePlotConfig, LinePlotConfig.from_df(df))

        # Reconstruct data structure using efficient grouping
        datasets: List[List[Series[KeyT, ValT]]] = []
        data_labels: List[Union[str, Dict[str, Any]]] = []
        sample_names: List[SampleName] = []
        sample_names_set: set = set()

        dataset_indices = sorted(df.select("dataset_idx").unique().to_series()) if not df.is_empty() else []

        for ds_idx in dataset_indices:
            ds_group = df.filter(pl.col("dataset_idx") == ds_idx) if not df.is_empty() else pl.DataFrame()

            data_label = ds_group.select("data_label").item(0, 0) if not ds_group.is_empty() else None
            data_labels.append(json.loads(data_label) if data_label else {})

            dataset: List[Series[KeyT, ValT]] = []

            if ds_group.is_empty():
                datasets.append(dataset)
                continue

            # Get unique sample names and sort them using natsort
            unique_samples_list = ds_group.select("sample").unique().to_series().to_list()
            sorted_samples = natsorted(unique_samples_list)

            # Build a lookup of sample -> rows using partition_by for efficiency
            # First, get all relevant columns as lists for faster access
            all_samples = ds_group.get_column("sample").to_list()
            all_x_vals = ds_group.get_column("x_val").to_list()
            all_y_vals = ds_group.get_column("y_val").to_list()
            all_x_types = ds_group.get_column("x_val_type").to_list()
            all_y_types = ds_group.get_column("y_val_type").to_list()
            all_series = ds_group.get_column("series").to_list()

            # Group data by sample name using a dictionary
            sample_data: Dict[str, List[int]] = {}
            for i, sample in enumerate(all_samples):
                if sample not in sample_data:
                    sample_data[sample] = []
                sample_data[sample].append(i)

            for sample_name in sorted_samples:
                row_indices = sample_data.get(sample_name, [])
                if not row_indices:
                    continue

                # Get series properties from first row
                first_idx = row_indices[0]
                series_dict = all_series[first_idx]

                # Extract x,y pairs
                pairs: List[Tuple[KeyT, ValT]] = []
                for idx in row_indices:
                    x_val = parse_value(all_x_vals[idx], all_x_types[idx])
                    y_val = parse_value(all_y_vals[idx], all_y_types[idx])
                    pairs.append((x_val, y_val))

                # Create Series object
                series: Series[KeyT, ValT] = Series(
                    name=str(sample_name),
                    pairs=pairs,
                    path_in_cfg=("lineplot", "data"),
                    **series_dict,
                )
                dataset.append(series)

                # Add sample name if not already in the set
                if sample_name not in sample_names_set:
                    sample_names_set.add(sample_name)
                    sample_names.append(SampleName(str(sample_name)))

            datasets.append(dataset)

        if any(d for d in data_labels if d):
            pconf.data_labels = data_labels
        return cls(
            anchor=anchor,
            plot_type=PlotType.LINE,
            data=datasets,
            pconfig=pconf,
            sample_names=sample_names,
            creation_date=cls.creation_date_from_df(df),
        )

    @staticmethod
    def create(
        data: Union[DatasetT[KeyT, ValT], Sequence[DatasetT[KeyT, ValT]]],
        pconfig: Union[Dict[str, Any], LinePlotConfig, None] = None,
    ) -> "LinePlotNormalizedInputData[KeyT, ValT]":
        pconf: LinePlotConfig = cast(LinePlotConfig, LinePlotConfig.from_pconfig_dict(pconfig))

        # Given one dataset - turn it into a list
        raw_dataset_list: List[DatasetT]
        if isinstance(data, Sequence):
            raw_dataset_list = list(data)
        else:
            raw_dataset_list = [data]
        del data

        # Normalise data labels
        if pconf.data_labels:
            if len(pconf.data_labels) != len(raw_dataset_list):
                raise ValueError(
                    f"Length of data_labels does not match the number of datasets. "
                    f"Please check your module code and ensure that the data_labels "
                    f"list is the same length as the data list: "
                    f"{len(pconf.data_labels)} != {len(raw_dataset_list)}. pconfig={pconf}"
                )

        sample_names = []

        datasets: List[List[Series[Any, Any]]] = []
        for ds_idx, raw_data_by_sample in enumerate(raw_dataset_list):
            list_of_series: List[Series[Any, Any]] = []
            for s in natsorted(raw_data_by_sample.keys()):
                if s not in sample_names:
                    sample_names.append(SampleName(s))
                x_to_y = raw_data_by_sample[s]
                if not isinstance(x_to_y, dict) and isinstance(x_to_y, Sequence):
                    if isinstance(x_to_y[0], tuple) or (isinstance(x_to_y[0], list) and len(x_to_y[0]) == 2):
                        x_to_y = dict(x_to_y)
                    else:
                        x_to_y = {i: y for i, y in enumerate(x_to_y)}
                dl = pconf.data_labels[ds_idx] if pconf.data_labels else None
                series: Series[Any, Any] = _make_series_dict(pconf, dl, s, x_to_y)
                if not series.pairs:
                    continue
                list_of_series.append(series)
            if list_of_series:
                datasets.append(list_of_series)

        # Return normalized data and config
        return LinePlotNormalizedInputData(
            anchor=plot_anchor(pconf),
            plot_type=PlotType.LINE,
            data=datasets,
            pconfig=pconf,
            sample_names=sample_names,
            creation_date=report.creation_date,
        )

    @classmethod
    def merge(
        cls,
        old_data: "LinePlotNormalizedInputData[KeyT, ValT]",
        new_data: "LinePlotNormalizedInputData[KeyT, ValT]",
    ) -> "LinePlotNormalizedInputData[KeyT, ValT]":
        """
        Merge normalized data from old run and new run, leveraging our tabular representation
        for more efficient and reliable merging.
        """
        logger.debug(f"LinePlot merge called for anchor {new_data.anchor}")
        new_df = new_data.to_df()
        if new_df.is_empty():
            logger.debug("LinePlot merge: new_df is empty, returning old_data")
            return old_data

        old_df = old_data.to_df()

        # If we have both old and new data, merge them
        merged_df = new_df
        if old_df is not None and not old_df.is_empty():
            logger.debug(
                f"LinePlot merge: both old and new data present, merging. Old: {old_df.shape}, New: {new_df.shape}"
            )
            # Get the list of samples that exist in both old and new data, for each dataset
            new_keys = new_df.select(["data_label", "sample"]).unique()
            logger.debug(f"LinePlot merge: new_keys shape: {new_keys.shape}")

            # Keep only the rows in old_df whose (data_label, sample) pair
            # does *not* appear in new_keys. An anti-join is the cleanest way to express this.
            old_df_filtered = old_df.join(
                new_keys,
                on=["data_label", "sample"],
                how="anti",  # anti-join = "rows in left not matched in right"
            )
            logger.debug(f"LinePlot merge: old_df_filtered shape: {old_df_filtered.shape}")

            # Combine the filtered old data with new data
            merged_df = pl.concat([old_df_filtered, new_df], how="diagonal")
            logger.debug(f"LinePlot merge: merged_df shape: {merged_df.shape}")
        else:
            logger.debug("LinePlot merge: no old data or old data is empty, using new data only")

        result = cls.from_df(merged_df, new_data.pconfig, new_data.anchor)
        logger.debug(
            f"LinePlot merge: created result with {len(result.data)} datasets, {len(result.sample_names)} samples"
        )
        return result


class LinePlot(Plot[Dataset[KeyT, ValT], LinePlotConfig], Generic[KeyT, ValT]):
    datasets: List[Dataset[KeyT, ValT]]
    sample_names: List[SampleName]

    def all_sample_names(self) -> List[SampleName]:
        return self.sample_names

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.xlab:
            result += f"X axis: {self.pconfig.xlab}\n"
        if self.pconfig.ylab:
            result += f"Y axis: {self.pconfig.ylab}\n"
        return result

    @staticmethod
    def create(
        lists_of_lines: List[List[Series[KeyT, ValT]]],
        pconfig: LinePlotConfig,
        anchor: Anchor,
        sample_names: List[SampleName],
    ) -> "LinePlot[KeyT, ValT]":
        lists_of_lines = [x for x in lists_of_lines if x]
        n_samples_per_dataset = [len(x) for x in lists_of_lines]

        axis_controlled_by_switches = pconfig.axis_controlled_by_switches or ["yaxis"]

        model: Plot[Dataset[KeyT, ValT], LinePlotConfig] = Plot.initialize(
            plot_type=PlotType.LINE,
            pconfig=pconfig,
            anchor=anchor,
            n_series_per_dataset=n_samples_per_dataset,
            axis_controlled_by_switches=list(axis_controlled_by_switches),
            default_tt_label="<br>%{x}: %{y}",
        )

        # Very large legend for automatically enabled flat plot mode is not very helpful
        max_n_samples = max(len(x) for x in lists_of_lines) if len(lists_of_lines) > 0 else 0
        if pconfig.showlegend is None and max_n_samples > 250:
            model.layout.showlegend = False

        model.datasets = [Dataset.create(d, lines, pconfig) for d, lines in zip(model.datasets, lists_of_lines)]

        # Make a tooltip always show on hover over any point on plot
        model.layout.hoverdistance = -1

        return LinePlot(**model.__dict__, sample_names=sample_names)

    @staticmethod
    def from_inputs(inputs: LinePlotNormalizedInputData[KeyT, ValT]) -> Union["LinePlot", str, None]:
        pconf = inputs.pconfig
        datasets = inputs.data
        sample_names = inputs.sample_names

        # Add extra annotation data series
        if pconf.extra_series:
            ess: Union[Series[Any, Any], List[Series[Any, Any]], List[List[Series[Any, Any]]]] = pconf.extra_series
            list_of_list_of_series: List[List[Series[Any, Any]]]
            if isinstance(ess, list):
                if isinstance(ess[0], list):
                    list_of_list_of_series = cast(List[List[Series[Any, Any]]], ess)
                else:
                    list_of_list_of_series = [cast(List[Series[Any, Any]], ess) for _ in datasets]
            else:
                list_of_list_of_series = [[ess] for _ in datasets]

            for i, list_of_raw_series in enumerate(list_of_list_of_series):
                assert isinstance(list_of_raw_series, list)
                for series in list_of_raw_series:
                    if i < len(datasets):
                        datasets[i].append(series)

        # Process categories
        for ds_idx, series_by_sample in enumerate(datasets):
            if pconf.categories and series_by_sample:
                if isinstance(pconf.categories, list):
                    categories = pconf.categories
                else:
                    categories = [pair[0] for pair in series_by_sample[0].pairs]
                for si, series in enumerate(series_by_sample):
                    if si != 0:
                        # If categories come in different order in different samples, reorder them
                        xs = [p[0] for p in series.pairs]
                        xs_set = set(xs)
                        xs_in_categories = [c for c in categories if c in xs_set]
                        categories_set = set(categories)
                        xs_not_in_categories = [x for x in xs if x not in categories_set]
                        xs = xs_in_categories + xs_not_in_categories
                        pairs = dict(series.pairs)
                        series.pairs = [(x, pairs[x]) for x in xs]

        inputs.save_to_parquet()

        scale = mqc_colour.mqc_colour_scale("plot_defaults")
        for _, series_by_sample in enumerate(datasets):
            for si, series in enumerate(series_by_sample):
                if not series.color:
                    series.color = scale.get_colour(si, lighten=1)

        plot = LinePlot.create(
            lists_of_lines=inputs.data,
            pconfig=inputs.pconfig,
            anchor=inputs.anchor,
            sample_names=sample_names,
        )
        return plot


def plot(
    data: Union[DatasetT[KeyT, ValT], Sequence[DatasetT[KeyT, ValT]]],
    pconfig: Union[Dict[str, Any], LinePlotConfig, None] = None,
) -> Union["LinePlot", str, None]:
    """
    Plot a line graph with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page

    Function effectively only returns a wrapper around parquet file path.
    """
    inputs: LinePlotNormalizedInputData[KeyT, ValT] = LinePlotNormalizedInputData.create(data, pconfig)
    inputs = LinePlotNormalizedInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None
    return LinePlot.from_inputs(inputs)


def remove_nones_and_empty_dicts(d: Mapping[Any, Any]) -> Dict[Any, Any]:
    """Remove None and empty dicts from a dict recursively."""
    return {k: remove_nones_and_empty_dicts(v) for k, v in d.items() if v is not None and v != {}}


def _make_series_dict(
    pconfig: LinePlotConfig,
    data_label: Union[Dict[str, Any], str, None],
    s: str,
    y_by_x: XToYDictT[KeyT, ValT],
) -> Series[KeyT, ValT]:
    pairs: List[Tuple[KeyT, ValT]] = []

    x_are_categories = pconfig.categories
    ymax = pconfig.ymax
    ymin = pconfig.ymin
    xmax = pconfig.xmax
    xmin = pconfig.xmin
    colors = pconfig.colors
    dash_styles = pconfig.dash_styles
    hovertemplates = pconfig.hovertemplates
    legend_groups = pconfig.legend_groups
    if data_label:
        if isinstance(data_label, dict):
            _x_are_categories = data_label.get("categories", x_are_categories)
            assert isinstance(_x_are_categories, bool)
            x_are_categories = _x_are_categories
            _ymax = data_label.get("ymax", ymax)
            _ymin = data_label.get("ymin", ymin)
            _xmax = data_label.get("xmax", xmax)
            _xmin = data_label.get("xmin", xmin)
            assert isinstance(_ymax, (int, float, type(None)))
            assert isinstance(_ymin, (int, float, type(None)))
            assert isinstance(_xmax, (int, float, type(None)))
            assert isinstance(_xmin, (int, float, type(None)))
            ymax = _ymax
            ymin = _ymin
            xmax = _xmax
            xmin = _xmin
            _colors = data_label.get("colors")
            if _colors and isinstance(_colors, dict):
                colors = {**colors, **cast(Dict[str, str], _colors)}
            _dash_styles = data_label.get("dash_styles")
            if _dash_styles and isinstance(_dash_styles, dict):
                dash_styles = {**dash_styles, **cast(Dict[str, str], _dash_styles)}
            _hovertemplates = data_label.get("hovertemplates")
            if _hovertemplates and isinstance(_hovertemplates, dict):
                hovertemplates = {**hovertemplates, **cast(Dict[str, str], _hovertemplates)}
            _legend_groups = data_label.get("legend_groups")
            if _legend_groups and isinstance(_legend_groups, dict):
                legend_groups = {**legend_groups, **cast(Dict[str, str], _legend_groups)}

    xs = [x for x in y_by_x.keys()]
    if not x_are_categories:
        xs = sorted(xs)

    # Discard > ymax or just hide?
    # If it never comes back into the plot, discard. If it goes above then comes back, just hide.
    discard_ymax = None
    discard_ymin = None
    for x in xs:
        if not x_are_categories:
            if xmax is not None and float(x) > float(xmax):
                continue
            if xmin is not None and float(x) < float(xmin):
                continue
        y = y_by_x[x]
        if y is not None:
            if ymax is not None:
                if float(y) > float(ymax):
                    discard_ymax = True
                elif discard_ymax is True:
                    discard_ymax = False
            if ymin is not None:
                if float(y) > float(ymin):
                    discard_ymin = True
                elif discard_ymin is True:
                    discard_ymin = False

    # Build the plot data structure
    for x in xs:
        if not x_are_categories and x is not None:
            if xmax is not None and float(x) > float(xmax):
                continue
            if xmin is not None and float(x) < float(xmin):
                continue

        y = y_by_x[x]
        if y is not None:
            if ymax is not None and float(y) > float(ymax) and discard_ymax is not False:
                continue
            if ymin is not None and float(y) < float(ymin) and discard_ymin is not False:
                continue
        pairs.append((x, y))

    # Smooth dataset if requested in config
    if pconfig.smooth_points is not None:
        pairs = smooth_array(pairs, pconfig.smooth_points)

    # Prepare extra trace parameters for hovertemplate and legendgroup
    extra_trace_params = {}
    hovertemplate = hovertemplates.get(s)
    if hovertemplate:
        extra_trace_params["hovertemplate"] = hovertemplate

    legendgroup = legend_groups.get(s)
    if legendgroup:
        extra_trace_params["legendgroup"] = legendgroup

    return Series(
        name=s,
        pairs=pairs,
        color=colors.get(s),
        dash=dash_styles.get(s),
        extra_trace_params=extra_trace_params,
        path_in_cfg=("lineplot", "pconfig", "pairs"),
    )


def smooth_line_data(data_by_sample: DatasetT[KeyT, ValT], numpoints: int) -> Dict[SampleName, Dict[KeyT, ValT]]:
    """
    Function to take an x-y dataset and use binning to smooth to a maximum number of datapoints.
    Each datapoint in a smoothed dataset corresponds to the first point in a bin.

    Examples to show the idea:

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=6
    we want to keep the first and the last element, thus excluding the last element from the binning:
    binsize = len([0 1 2 3 4 5 6 7 8]))/(numpoints-1) = 9/5 = 1.8
    taking points in indices rounded from multiples of 1.8: [0, 1.8, 3.6, 5.4, 7.2, 9],
    ...which evaluates to first_element_in_bin_indices=[0, 2, 4, 5, 7, 9]
    picking up the elements: [0 _ 2 _ 4 5 _ 7 _ 9]

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=9
    binsize = 9/8 = 1.125
    indices: [0.0, 1.125, 2.25, 3.375, 4.5, 5.625, 6.75, 7.875, 9] -> [0, 1, 2, 3, 5, 6, 7, 8, 9]
    picking up the elements: [0 1 2 3 _ 5 6 7 8 9]

    d=[0 1 2 3 4 5 6 7 8 9], numpoints=3
    binsize = len(d)/numpoints = 9/2 = 4.5
    indices: [0.0, 4.5, 9] -> [0, 5, 9]
    picking up the elements: [0 _ _ _ _ 5 _ _ _ 9]
    """
    smoothed_data: Dict[SampleName, Dict[KeyT, ValT]] = dict()
    for s_name, d in data_by_sample.items():
        smoothed_data[SampleName(s_name)] = dict(smooth_array(list(d.items()), numpoints))

    return smoothed_data


T = TypeVar("T")


def smooth_array(items: List[T], numpoints: int) -> List[T]:
    """
    Function to take an array and use binning to smooth to a maximum number of datapoints.
    Each datapoint in a smoothed dataset corresponds to the first point in a bin.
    """
    # Check that we need to smooth this data
    if len(items) <= numpoints or len(items) == 0:
        return items

    result: List[T] = []
    binsize = (len(items) - 1) / (numpoints - 1)
    first_element_indices = {round(binsize * i) for i in range(numpoints)}
    for i, y in enumerate(items):
        if i in first_element_indices:
            result.append(y)
    return result
