"""MultiQC functions to plot a box plot"""

import copy
import json
import logging
from typing import Any, Dict, List, Mapping, Optional, OrderedDict, Tuple, Union, cast

import plotly.graph_objects as go  # type: ignore
import polars as pl
from natsort import natsorted

from multiqc import config, report
from multiqc.plots.plot import BaseDataset, NormalizedPlotInputData, PConfig, Plot, PlotType, plot_anchor
from multiqc.plots.utils import determine_barplot_height
from multiqc.types import Anchor, SampleName

logger = logging.getLogger(__name__)


class BoxPlotConfig(PConfig):
    sort_samples: bool = True
    sort_by_median: bool = False  # Sort samples by their median values instead of alphabetically
    sort_switch_sorted_active: bool = True  # Whether the sorted view is initially active
    # Showing separate points:
    # False - do not show at all
    # "all" - show all data points
    # "outliers" - show only outliers
    # None - use config.boxplot_boxpoints; if not set, determine based on config.box_min_threshold_no_points and config.box_min_threshold_outliers
    boxpoints: Union[bool, str, None] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("boxplot",), **data)


# Type of single box (matching one sample) - can be raw data or statistics
BoxT = Union[List[Union[int, float]], Dict[str, Union[int, float]]]
# Type for statistics dict
BoxStatsT = Dict[str, Union[int, float]]


class Dataset(BaseDataset):
    data: List[BoxT]
    samples: List[str]
    data_sorted: Optional[List[BoxT]] = None  # Sorted version of data
    samples_sorted: Optional[List[str]] = None  # Sorted version of samples
    is_stats_data: bool = False  # True if data contains pre-calculated statistics

    def sample_names(self) -> List[SampleName]:
        return [SampleName(sample) for sample in self.samples]

    @staticmethod
    def create(
        dataset: BaseDataset,
        data_by_sample: Mapping[str, BoxT],
        pconfig: Optional[BoxPlotConfig] = None,
    ) -> "Dataset":
        # Detect if we have statistics data or raw data
        is_stats_data = False
        if data_by_sample:
            first_sample_data = next(iter(data_by_sample.values()))
            if isinstance(first_sample_data, dict):
                is_stats_data = True
        # Store original order (reversed for box plot display)
        original_data = list(data_by_sample.values())
        original_samples = list(data_by_sample.keys())
        original_data = list(reversed(original_data))
        original_samples = list(reversed(original_samples))

        # Store sorted version if sort_by_median is enabled
        data_sorted = None
        samples_sorted = None
        main_data = original_data
        main_samples = original_samples

        if pconfig and pconfig.sort_by_median:
            # Calculate median for each sample and sort by it
            median_values = {}
            for sample, values in data_by_sample.items():
                if values:
                    if is_stats_data:
                        # Use pre-calculated median from statistics
                        stats_dict = cast(BoxStatsT, values)
                        median_values[sample] = stats_dict.get("median", 0)
                    else:
                        # Calculate median from raw data
                        raw_values = cast(List[Union[int, float]], values)
                        sorted_values = sorted(raw_values)
                        n = len(sorted_values)
                        median = (
                            sorted_values[n // 2]
                            if n % 2 == 1
                            else (sorted_values[n // 2 - 1] + sorted_values[n // 2]) / 2
                        )
                        median_values[sample] = median
                else:
                    median_values[sample] = 0  # Handle empty data

            # Sort samples by median (descending order for better visual hierarchy)
            sorted_sample_names = sorted(median_values.keys(), key=lambda x: median_values[x], reverse=True)
            sorted_data_by_sample = {sample: data_by_sample[sample] for sample in sorted_sample_names}

            # Store sorted versions (also reversed for box plot display)
            data_sorted = list(reversed(list(sorted_data_by_sample.values())))
            samples_sorted = list(reversed(list(sorted_data_by_sample.keys())))

            # If starting with sorted view active, keep sorted data in data_sorted
            # JavaScript will choose based on sortSwitchSortedActive flag
            if pconfig.sort_switch_sorted_active:
                # Keep the original data as main, sorted data as data_sorted
                # JavaScript will use data_sorted when sortSwitchSortedActive is true
                main_data = original_data
                main_samples = original_samples

        dataset = Dataset(
            **dataset.__dict__,
            data=main_data,
            samples=main_samples,
            data_sorted=data_sorted,
            samples_sorted=samples_sorted,
            is_stats_data=is_stats_data,
        )

        # Determine boxpoints based on PConfig first, then global config, then dynamic logic
        boxpoints: Union[bool, str] = "outliers"

        if is_stats_data:
            # For statistics data, we can't show individual points
            boxpoints = False
        elif pconfig and pconfig.boxpoints is not None:
            # Use explicit PConfig boxpoints setting
            boxpoints = pconfig.boxpoints
        elif config.boxplot_boxpoints is not None:
            # Use global config setting
            boxpoints = config.boxplot_boxpoints
        else:
            # Fall back to dynamic logic based on sample count, similar to violin plot
            n_samples = len(dataset.samples)
            show_points = n_samples <= config.box_min_threshold_no_points
            show_only_outliers = n_samples > config.box_min_threshold_outliers

            # Set boxpoints based on dynamic logic
            if not show_points:
                boxpoints = False  # Show no points
            elif not show_only_outliers:
                boxpoints = "all"  # Show all points
            else:
                boxpoints = "outliers"  # Show only outliers

        dataset.trace_params.update(
            boxpoints=boxpoints,
            jitter=0.5,
            orientation="h",
            marker=dict(
                color="#4899e8",  # just use blue to indicate interactivity
            ),
            # to remove the redundant sample name before "median" in the unified hover box
            hoverinfo="x",
        )
        dataset.layout["yaxis"]["title"] = None
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log: bool = False,  # noqa: ARG002
        is_pct: bool = False,  # noqa: ARG002
        **kwargs,  # noqa: ARG002
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)

        for sname, values in zip(self.samples, self.data):
            params = copy.deepcopy(self.trace_params)

            if self.is_stats_data:
                # Use statistics to create box plot
                stats_dict = cast(BoxStatsT, values)
                fig.add_trace(
                    go.Box(
                        q1=[stats_dict.get("q1", 0)],
                        median=[stats_dict.get("median", 0)],
                        q3=[stats_dict.get("q3", 0)],
                        lowerfence=[stats_dict.get("min", 0)],
                        upperfence=[stats_dict.get("max", 0)],
                        mean=[stats_dict.get("mean", stats_dict.get("median", 0))],
                        name=sname,
                        **params,
                    ),
                )
            else:
                # Use raw data points
                raw_values = cast(List[Union[int, float]], values)
                fig.add_trace(
                    go.Box(
                        x=raw_values,
                        name=sname,
                        **params,
                    ),
                )
        return fig

    def save_data_file(self) -> None:
        vals_by_sample: Dict[str, BoxT] = {}
        for sample, values in zip(self.samples, self.data):
            if self.is_stats_data:
                # For statistics data, save the statistics dict
                vals_by_sample[sample] = values
            else:
                # For raw data, save the raw values
                vals_by_sample[sample] = values
        report.write_data_file(vals_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:  # type: ignore[override]
        """Format dataset as a markdown table with basic statistics"""
        prompt = "|Sample|Min|Q1|Median|Q3|Max|Mean|\n"
        prompt += "|---|---|---|---|---|---|---|\n"

        suffix = ""
        if self.layout["xaxis"]["ticksuffix"]:
            suffix = " " + self.layout["xaxis"]["ticksuffix"]

        for sample, values in zip(self.samples, self.data):
            # Skip samples with no data
            if (self.is_stats_data and not values) or (not self.is_stats_data and len(values) == 0):
                continue

            # Use pseudonym if available, otherwise use original sample name
            pseudonym = report.anonymize_sample_name(sample)

            if self.is_stats_data:
                # Use pre-calculated statistics
                stats_dict = cast(BoxStatsT, values)
                min_val = stats_dict.get("min", 0)
                max_val = stats_dict.get("max", 0)
                median = stats_dict.get("median", 0)
                q1 = stats_dict.get("q1", min_val)
                q3 = stats_dict.get("q3", max_val)
                mean = stats_dict.get("mean", median)
            else:
                # Calculate statistics from raw data
                raw_values = cast(List[Union[int, float]], values)
                sorted_vals = sorted(raw_values)
                n = len(sorted_vals)

                min_val = sorted_vals[0]
                max_val = sorted_vals[-1]
                median = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
                q1 = sorted_vals[n // 4] if n >= 4 else sorted_vals[0]
                q3 = sorted_vals[3 * n // 4] if n >= 4 else sorted_vals[-1]
                mean = sum(raw_values) / len(raw_values)

            prompt += (
                f"|{pseudonym}|"
                f"{self.fmt_value_for_llm(min_val)}{suffix}|"
                f"{self.fmt_value_for_llm(q1)}{suffix}|"
                f"{self.fmt_value_for_llm(median)}{suffix}|"
                f"{self.fmt_value_for_llm(q3)}{suffix}|"
                f"{self.fmt_value_for_llm(max_val)}{suffix}|"
                f"{self.fmt_value_for_llm(mean)}{suffix}|\n"
            )
        return prompt


class BoxPlotInputData(NormalizedPlotInputData):
    list_of_data_by_sample: List[Mapping[str, BoxT]]
    pconfig: BoxPlotConfig

    def is_empty(self) -> bool:
        return len(self.list_of_data_by_sample) == 0 or all(len(ds) == 0 for ds in self.list_of_data_by_sample)

    def to_df(self) -> pl.DataFrame:
        """
        Save plot data to a parquet file using a tabular representation that's
        optimized for cross-run analysis.
        """
        if self.is_empty():
            return pl.DataFrame()

        records = []

        for ds_idx, dataset in enumerate(self.list_of_data_by_sample):
            # Process each sample in the dataset
            for sample_name, values in dataset.items():
                # Store each value point as a separate record
                for value_idx, value in enumerate(values):
                    record = {
                        "anchor": self.anchor,
                        "dataset_idx": ds_idx,
                        "data_label": json.dumps(self.pconfig.data_labels[ds_idx]) if self.pconfig.data_labels else "",
                        "sample": str(sample_name),
                        "value_idx": value_idx,
                        "value": value,
                    }

                    # Add plot configuration metadata
                    if hasattr(self.pconfig, "title") and self.pconfig.title:
                        record["plot_title"] = self.pconfig.title
                    if hasattr(self.pconfig, "xlab") and self.pconfig.xlab:
                        record["x_label"] = self.pconfig.xlab
                    if hasattr(self.pconfig, "ylab") and self.pconfig.ylab:
                        record["y_label"] = self.pconfig.ylab

                    records.append(record)

        df = pl.DataFrame(records, schema_overrides={"data_label": pl.Utf8})
        return self.finalize_df(df)

    @classmethod
    def from_df(cls, df: pl.DataFrame, pconfig: Union[Dict, BoxPlotConfig], anchor: Anchor) -> "BoxPlotInputData":
        """
        Load plot data from a DataFrame.
        """
        if df.is_empty():
            pconf = (
                pconfig
                if isinstance(pconfig, BoxPlotConfig)
                else cast(BoxPlotConfig, BoxPlotConfig.from_pconfig_dict(pconfig))
            )
            return cls(
                anchor=anchor,
                pconfig=pconf,
                list_of_data_by_sample=[],
                plot_type=PlotType.BOX,
                creation_date=cls.creation_date_from_df(df),
            )
        pconf = cast(BoxPlotConfig, BoxPlotConfig.from_df(df))

        # Group by dataset_idx to rebuild data structure
        list_of_data_by_sample: List[Mapping[str, BoxT]] = []
        data_labels = []

        max_dataset_idx = df.select(pl.col("dataset_idx").max()).item() if not df.is_empty() else 0
        for ds_idx in range(int(max_dataset_idx) + 1):
            ds_group = df.filter(pl.col("dataset_idx") == ds_idx) if not df.is_empty() else pl.DataFrame()

            data_label = ds_group.select("data_label").item(0) if not ds_group.is_empty() else None
            data_labels.append(json.loads(data_label) if data_label else {})

            # Process each sample in this dataset
            dataset = {}

            unique_samples: pl.Series = (
                ds_group.select("sample").unique().to_series() if not ds_group.is_empty() else pl.Series([])
            )
            # Group by sample_name within each dataset
            for sample_name in unique_samples:
                # Get data for this sample
                sample_group = ds_group.filter(pl.col("sample") == sample_name)

                # Collect all values for this sample
                sample_values = []

                # Sort by value_idx if available
                if "value_idx" in sample_group.columns:
                    sample_group = sample_group.sort("value_idx")

                for row in sample_group.iter_rows(named=True):
                    value = row["value"]
                    sample_values.append(value)

                # Add sample data if not empty
                if sample_values:
                    dataset[str(sample_name)] = sample_values

            list_of_data_by_sample.append(cast(Mapping[str, BoxT], dataset))

        if any(d for d in data_labels if d):
            pconf.data_labels = data_labels
        return cls(
            anchor=anchor,
            pconfig=pconf,
            list_of_data_by_sample=list_of_data_by_sample,
            plot_type=PlotType.BOX,
            creation_date=cls.creation_date_from_df(df),
        )

    @classmethod
    def merge(cls, old_data: "BoxPlotInputData", new_data: "BoxPlotInputData") -> "BoxPlotInputData":
        """
        Merge normalized data from old run and new run, matching by data labels when available
        """
        # Create dataframe for new data
        new_df = new_data.to_df()
        if new_df.is_empty():
            return old_data

        old_df = old_data.to_df()

        # If we have both old and new data, merge them
        merged_df = new_df
        if old_df is not None and not old_df.is_empty():
            # Combine the dataframes, keeping all rows
            merged_df = pl.concat([old_df, new_df], how="diagonal")

            # For duplicates (same sample, dataset, value_idx), keep the latest version
            # Sort by timestamp (newest last)
            merged_df = merged_df.sort("creation_date")

            # Group by the key identifiers and keep the last entry (newest)
            dedupe_columns = ["dataset_idx", "sample", "value_idx"]
            merged_df = (
                merged_df.group_by(dedupe_columns).agg(pl.all().exclude(dedupe_columns).last()).sort("creation_date")
            )

        # Create a new BoxPlotInputData from the merged DataFrame
        return cls.from_df(merged_df, new_data.pconfig, new_data.anchor)

    @staticmethod
    def create(
        list_of_data_by_sample: Union[Mapping[str, BoxT], List[Mapping[str, BoxT]]],
        pconfig: Union[Dict[str, Any], BoxPlotConfig, None] = None,
    ) -> "BoxPlotInputData":
        pconf: BoxPlotConfig = cast(BoxPlotConfig, BoxPlotConfig.from_pconfig_dict(pconfig))

        # Given one dataset - turn it into a list
        if not isinstance(list_of_data_by_sample, list):
            list_of_data_by_sample = [list_of_data_by_sample]

        for i in range(len(list_of_data_by_sample)):
            if isinstance(list_of_data_by_sample[0], OrderedDict):
                # Legacy: users assumed that passing an OrderedDict indicates that we
                # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
                pass
            elif pconf.sort_samples:
                samples = natsorted(list(list_of_data_by_sample[0].keys()))
                list_of_data_by_sample[i] = {s: list_of_data_by_sample[i][s] for s in samples}

        return BoxPlotInputData(
            anchor=plot_anchor(pconf),
            list_of_data_by_sample=list_of_data_by_sample,
            pconfig=pconf,
            plot_type=PlotType.BOX,
            creation_date=report.creation_date,
        )


class BoxPlot(Plot[Dataset, BoxPlotConfig]):
    datasets: List[Dataset]
    sort_switch_sorted_active: bool = False

    def sample_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        for ds in self.datasets:
            names.extend(ds.sample_names())
        return names

    @staticmethod
    def create(
        list_of_data_by_sample: List[Mapping[str, BoxT]],
        pconfig: BoxPlotConfig,
        anchor: Anchor,
    ) -> "BoxPlot":
        model: Plot[Dataset, BoxPlotConfig] = Plot.initialize(
            plot_type=PlotType.BOX,
            pconfig=pconfig,
            anchor=anchor,
            n_series_per_dataset=[len(x) for x in list_of_data_by_sample],
        )

        model.datasets = [
            Dataset.create(ds, data_by_sample, pconfig)
            for ds, data_by_sample in zip(model.datasets, list_of_data_by_sample)
        ]

        max_n_samples = max(len(x) for x in list_of_data_by_sample) if list_of_data_by_sample else 0
        height: int = determine_barplot_height(max_n_samples)

        model.layout.update(
            height=height,
            showlegend=False,
            boxgroupgap=0.1,
            boxgap=0.2,
            colorway=[],  # no need to color code
            yaxis=dict(
                automargin=True,  # to make sure there is enough space for ticks labels
                categoryorder="trace",  # keep sample order
                hoverformat=getattr(model.layout.xaxis, "hoverformat", None),
                ticksuffix=getattr(model.layout.xaxis, "ticksuffix", None),
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                title=dict(text=getattr(getattr(model.layout.yaxis, "title", None), "text", None)),
                hoverformat=getattr(model.layout.yaxis, "hoverformat", None),
                ticksuffix=getattr(model.layout.yaxis, "ticksuffix", None),
            ),
            hovermode="y",
            hoverlabel=dict(
                bgcolor="white",
                font=dict(color="rgba(60,60,60,1)"),
            ),
        )
        return BoxPlot(**model.__dict__, sort_switch_sorted_active=pconfig.sort_switch_sorted_active)

    def buttons(self, flat: bool, module_anchor: Anchor, section_anchor: Anchor) -> List[str]:
        """
        Box plot-specific controls, only for the interactive version.
        """
        buttons = super().buttons(flat=flat, module_anchor=module_anchor, section_anchor=section_anchor)
        if any(ds.data_sorted for ds in self.datasets):
            buttons.append(
                f"""
                <div class="btn-group" role="group">
                    <button
                        type="button"
                        class="btn btn-default btn-sm {"" if self.pconfig.sort_switch_sorted_active else "active"}"
                        data-action="unsorted"
                        data-plot-anchor="{self.anchor}"
                    >
                        Sorted by label
                    </button>
                    <button
                        type="button"
                        class="btn btn-default btn-sm {"active" if self.pconfig.sort_switch_sorted_active else ""}"
                        data-action="sorted_by_median"
                        data-plot-anchor="{self.anchor}"
                    >
                        Sorted by median
                    </button>
                </div>
                """
            )
        return buttons

    @staticmethod
    def from_inputs(inputs: BoxPlotInputData) -> Union["BoxPlot", str, None]:
        plot = BoxPlot.create(
            list_of_data_by_sample=inputs.list_of_data_by_sample,
            pconfig=inputs.pconfig,
            anchor=inputs.anchor,
        )
        inputs.save_to_parquet()
        return plot


def plot(
    list_of_data_by_sample: Union[Mapping[str, BoxT], List[Mapping[str, BoxT]]],
    pconfig: Union[Dict[str, Any], BoxPlotConfig, None] = None,
) -> Union["BoxPlot", str, None]:
    """
    Plot a box plot. Supports two input formats:

    1. Raw data points (traditional):
       {'sample1': [1.2, 3.4, 2.1, ...], 'sample2': [2.3, 4.1, ...]}

    2. Pre-calculated statistics (memory efficient for large datasets):
       {'sample1': {'min': 1.0, 'q1': 2.0, 'median': 3.0, 'q3': 4.0, 'max': 5.0, 'mean': 3.2},
        'sample2': {'min': 1.5, 'q1': 2.2, 'median': 3.1, 'q3': 4.1, 'max': 5.2, 'mean': 3.3}}

    The statistics format dramatically reduces memory usage and file sizes for large datasets
    (e.g., cell-level measurements) while producing identical visual output.

    Required statistics keys: min, q1, median, q3, max
    Optional statistics keys: mean, count, std
    """
    inputs: BoxPlotInputData = BoxPlotInputData.create(list_of_data_by_sample, pconfig)
    inputs = BoxPlotInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    return BoxPlot.from_inputs(inputs)
