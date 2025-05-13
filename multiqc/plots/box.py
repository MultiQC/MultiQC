"""MultiQC functions to plot a box plot"""

import copy
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, OrderedDict, Tuple, Union, cast

import pandas as pd
import plotly.graph_objects as go  # type: ignore

from multiqc import config, report
from multiqc.core.plot_data_store import append_to_parquet
from multiqc.plots.plot import BaseDataset, NormalizedPlotInputData, PConfig, Plot, PlotType, plot_anchor
from multiqc.plots.utils import determine_barplot_height
from multiqc.types import Anchor, SampleName

logger = logging.getLogger(__name__)


class BoxPlotConfig(PConfig):
    sort_samples: bool = True

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("boxplot",), **data)


# Type of single box (matching one sample)
BoxT = List[Union[int, float]]


class Dataset(BaseDataset):
    data: List[BoxT]
    samples: List[str]

    def sample_names(self) -> List[SampleName]:
        return [SampleName(sample) for sample in self.samples]

    @staticmethod
    def create(
        dataset: BaseDataset,
        data_by_sample: Dict[str, BoxT],
    ) -> "Dataset":
        dataset = Dataset(
            **dataset.__dict__,
            data=list(data_by_sample.values()),
            samples=list(data_by_sample.keys()),
        )
        # Need to reverse samples as the box plot will show them reversed
        dataset.samples = list(reversed(dataset.samples))
        dataset.data = list(reversed(dataset.data))

        dataset.trace_params.update(
            boxpoints="outliers",
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
        is_log: bool = False,
        is_pct: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)

        for sname, values in zip(self.samples, self.data):
            params = copy.deepcopy(self.trace_params)
            fig.add_trace(
                go.Box(
                    x=values,
                    name=sname,
                    **params,
                ),
            )
        return fig

    def save_data_file(self) -> None:
        vals_by_sample: Dict[str, BoxT] = {}
        for sample, values in zip(self.samples, self.data):
            vals_by_sample[sample] = values
        report.write_data_file(vals_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        """Format dataset as a markdown table with basic statistics"""
        prompt = "|Sample|Min|Q1|Median|Q3|Max|Mean|\n"
        prompt += "|---|---|---|---|---|---|---|\n"

        suffix = ""
        if self.layout["xaxis"]["ticksuffix"]:
            suffix = " " + self.layout["xaxis"]["ticksuffix"]

        for sample, values in zip(self.samples, self.data):
            # Skip samples with no data
            if len(values) == 0:
                continue

            # Use pseudonym if available, otherwise use original sample name
            pseudonym = report.anonymize_sample_name(sample)

            # Calculate statistics
            sorted_vals = sorted(values)
            n = len(sorted_vals)

            min_val = sorted_vals[0]
            max_val = sorted_vals[-1]
            median = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
            q1 = sorted_vals[n // 4] if n >= 4 else sorted_vals[0]
            q3 = sorted_vals[3 * n // 4] if n >= 4 else sorted_vals[-1]
            mean = sum(values) / len(values)

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
    list_of_data_by_sample: List[Dict[str, BoxT]]
    pconfig: BoxPlotConfig

    def is_empty(self) -> bool:
        return len(self.list_of_data_by_sample) == 0 or all(len(ds) == 0 for ds in self.list_of_data_by_sample)

    def to_df(self) -> pd.DataFrame:
        """
        Save plot data to a parquet file using a tabular representation that's
        optimized for cross-run analysis.
        """
        if self.is_empty():
            return pd.DataFrame()

        records = []

        for ds_idx, dataset in enumerate(self.list_of_data_by_sample):
            # Process each sample in the dataset
            for sample_name, values in dataset.items():
                # Store each value point as a separate record
                for value_idx, value in enumerate(values):
                    record = {
                        "anchor": self.anchor,
                        "dataset_idx": ds_idx,
                        "data_label": json.dumps(self.pconfig.data_labels[ds_idx])
                        if self.pconfig.data_labels
                        else None,
                        "sample_name": str(sample_name),
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

        df = pd.DataFrame(records)
        self.finalize_df(df)
        return df

    @classmethod
    def from_df(cls, df: pd.DataFrame, pconfig: Union[Dict, BoxPlotConfig], anchor: Anchor) -> "BoxPlotInputData":
        """
        Load plot data from a DataFrame.
        """
        if df.empty:
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
        list_of_data_by_sample: List[Dict[str, BoxT]] = []
        data_labels = []

        max_dataset_idx = df["dataset_idx"].max() if not df.empty else 0
        for ds_idx in range(int(max_dataset_idx) + 1):
            ds_group = df[df["dataset_idx"] == ds_idx] if not df.empty else pd.DataFrame()
            data_label = ds_group["data_label"].iloc[0]
            data_labels.append(json.loads(data_label) if data_label else {})

            # Process each sample in this dataset
            dataset = {}

            unique_samples = ds_group["sample_name"].unique()
            # Group by sample_name within each dataset
            for sample_name in unique_samples:
                # Get data for this sample
                sample_group = ds_group[ds_group["sample_name"] == sample_name]

                # Collect all values for this sample
                sample_values = []

                # Sort by value_idx if available
                if "value_idx" in sample_group.columns:
                    sample_group = sample_group.sort_values("value_idx")

                for _, row in sample_group.iterrows():
                    value = row["value"]
                    sample_values.append(value)

                # Add sample data if not empty
                if sample_values:
                    dataset[str(sample_name)] = sample_values

            list_of_data_by_sample.append(dataset)

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
        if new_df.empty:
            return old_data

        old_df = old_data.to_df()

        # If we have both old and new data, merge them
        merged_df = new_df
        if old_df is not None and not old_df.empty:
            # Combine the dataframes, keeping all rows
            merged_df = pd.concat([old_df, new_df], ignore_index=True)

            # For duplicates (same sample, dataset, value_idx), keep the latest version
            # Sort by timestamp (newest last)
            merged_df.sort_values("creation_date", inplace=True)

            # Group by the key identifiers and keep the last entry (newest)
            dedupe_columns = ["dataset_idx", "sample_name", "value_idx"]
            merged_df = merged_df.drop_duplicates(subset=dedupe_columns, keep="last")

        # Save the merged data for future reference
        append_to_parquet(merged_df)

        # Create a new BoxPlotInputData from the merged DataFrame
        return cls.from_df(merged_df, new_data.pconfig, new_data.anchor)

    @staticmethod
    def create(
        list_of_data_by_sample: Union[Dict[str, BoxT], List[Dict[str, BoxT]]],
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
                samples = sorted(list(list_of_data_by_sample[0].keys()))
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

    def sample_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        for ds in self.datasets:
            names.extend(ds.sample_names())
        return names

    @staticmethod
    def create(
        list_of_data_by_sample: List[Dict[str, BoxT]],
        pconfig: BoxPlotConfig,
        anchor: Anchor,
    ) -> "BoxPlot":
        model: Plot[Dataset, BoxPlotConfig] = Plot.initialize(
            plot_type=PlotType.BOX,
            pconfig=pconfig,
            anchor=anchor,
            n_samples_per_dataset=[len(x) for x in list_of_data_by_sample],
        )

        model.datasets = [
            Dataset.create(ds, data_by_sample) for ds, data_by_sample in zip(model.datasets, list_of_data_by_sample)
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
                hoverformat=model.layout.xaxis.hoverformat,
                ticksuffix=model.layout.xaxis.ticksuffix,
                # Prevent JavaScript from automatically parsing categorical values as numbers:
                type="category",
            ),
            xaxis=dict(
                title=dict(text=model.layout.yaxis.title.text),
                hoverformat=model.layout.yaxis.hoverformat,
                ticksuffix=model.layout.yaxis.ticksuffix,
            ),
            hovermode="y",
            hoverlabel=dict(
                bgcolor="rgba(255, 255, 255, 0.8)",
                font=dict(color="black"),
            ),
        )
        return BoxPlot(**model.__dict__)

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
    list_of_data_by_sample: Union[Dict[str, BoxT], List[Dict[str, BoxT]]],
    pconfig: Union[Dict[str, Any], BoxPlotConfig, None] = None,
) -> Union["BoxPlot", str, None]:
    """
    Plot a box plot. Expects either:
    - a dict mapping sample names to data point lists or dicts,
    - a dict mapping sample names to a dict of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    """
    inputs: BoxPlotInputData = BoxPlotInputData.create(list_of_data_by_sample, pconfig)
    inputs = BoxPlotInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    return BoxPlot.from_inputs(inputs)
