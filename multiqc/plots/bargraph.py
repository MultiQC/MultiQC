"""MultiQC functions to plot a bargraph"""

import copy
import logging
import math
from collections import OrderedDict, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Literal, Mapping, NewType, Optional, Sequence, Tuple, TypedDict, Union, cast

import pandas as pd
import plotly.graph_objects as go  # type: ignore
import spectra  # type: ignore
from natsort import natsorted
from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.core import tmp_dir
from multiqc.core.exceptions import RunError
from multiqc.plots.plot import (
    BaseDataset,
    NormalizedPlotInputData,
    PConfig,
    Plot,
    PlotType,
    plot_anchor,
    split_long_string,
)
from multiqc.plots.utils import determine_barplot_height
from multiqc.types import Anchor, SampleName
from multiqc.utils import mqc_colour
from multiqc.validation import ValidatedConfig

logger = logging.getLogger(__name__)

SampleNameT = Union[SampleName, str]
CatName = NewType("CatName", str)
CatNameT = Union[CatName, str]
InputDatasetT = Union[Mapping[SampleName, Mapping[CatName, Any]], Mapping[str, Mapping[str, Any]]]


class CatConf(ValidatedConfig):
    name: str
    color: Optional[str] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("cats",), **data)


# Either a list of strings, or a cat conf - a mapping from category names to their properties dicts or objects
InputCategoriesT = Union[
    Sequence[str],
    Mapping[CatName, Mapping[str, Any]],
    Sequence[CatName],
    Mapping[str, Mapping[str, Any]],
    Sequence[str],
]


class BarPlotConfig(PConfig):
    stacking: Union[Literal["group", "overlay", "relative", "normal"], None] = "relative"
    hide_empty: Optional[bool] = Field(None, deprecated="hide_zero_cats")
    hide_zero_cats: bool = True
    sort_samples: bool = True
    use_legend: Optional[bool] = None
    suffix: Optional[str] = None
    lab_format: Optional[str] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        if "suffix" in data:
            data["ysuffix"] = data["suffix"]
            del data["suffix"]
        if "lab_format" in data:
            data["ylab_format"] = data["lab_format"]
            del data["lab_format"]

        super().__init__(path_in_cfg=path_in_cfg or ("barplot",), **data)


class CatDataDict(TypedDict):
    """
    The way data is prepared and serialized for plotting
    """

    name: str
    color: str
    data: List[float]
    data_pct: List[float]


DatasetT = Dict[SampleName, Dict[CatName, Union[int, float]]]


class BarPlotInputData(NormalizedPlotInputData):
    pconfig: BarPlotConfig
    data: List[DatasetT]
    cats: List[Dict[CatName, CatConf]]

    def is_empty(self) -> bool:
        return len(self.data) == 0 or all(len(ds) == 0 for ds in self.data)

    @staticmethod
    def create(
        data: Union[InputDatasetT, Sequence[InputDatasetT]],
        cats: Optional[Union[InputCategoriesT, Sequence[InputCategoriesT]]] = None,
        pconfig: Optional[Union[Dict[str, Any], BarPlotConfig]] = None,
    ) -> "BarPlotInputData":
        """
        We want to be permissive with user input, e.g. allow one dataset or a list of datasets,
        allow optional categories, categories as strings or as dicts or as objects. We want
        to normalize the input data before we save it to intermediate format and plot.
        """
        pconf = cast(BarPlotConfig, BarPlotConfig.from_pconfig_dict(pconfig))

        # Given one dataset - turn it into a list
        raw_datasets: List[DatasetT]
        if isinstance(data, Sequence):
            raw_datasets = cast(List[DatasetT], data)
        else:
            raw_datasets = [cast(DatasetT, data)]
        del data

        # Make list of cats from different inputs
        raw_cats_per_ds: List[InputCategoriesT]
        if not cats:
            # Not supplied, generate default categories
            raw_cats_per_ds = []
            for val_by_cat_by_sample in raw_datasets:
                ds_cats: List[CatName] = []
                for sample_name, val_by_cat in val_by_cat_by_sample.items():
                    for _cat_name in val_by_cat.keys():
                        if _cat_name not in raw_cats_per_ds:
                            ds_cats.append(CatName(_cat_name))
                raw_cats_per_ds.append(ds_cats)
        elif isinstance(cats, List) and isinstance(cats[0], str):
            # ["Cat1", "Cat2"] - list of strings for one dataset
            raw_cats_per_ds = [[CatName(cat_name) for cat_name in cast(List[str], cats)]]
        elif isinstance(cats, Sequence):
            # [["Cat1", "Cat2"], {"Cat3": {}, "Cat4": {}}] - list of lists or dicts for multiple datasets
            raw_cats_per_ds = [ds_cats for ds_cats in cast(List[Dict], cats)]
        else:
            raw_cats_per_ds = [cats]

        if len(raw_datasets) > 1 and len(raw_cats_per_ds) == 1:
            raw_cats_per_ds = raw_cats_per_ds * len(raw_datasets)
        elif len(raw_datasets) == 0 and len(raw_cats_per_ds) == 1:
            raw_cats_per_ds = []
        elif len(raw_datasets) != len(raw_cats_per_ds):
            raise RunError(
                f"Bar graph: number of dataset and category lists must match, got {len(raw_datasets)} "
                f"datasets and {len(raw_cats_per_ds)} category lists: {raw_cats_per_ds}"
            )

        # Parse the categories into pydantic objects
        categories_per_ds: List[Dict[CatName, CatConf]] = []
        for raw_ds_cats in raw_cats_per_ds:
            ds_categories: Dict[CatName, CatConf] = dict()
            if isinstance(raw_ds_cats, list):
                for cat_name in raw_ds_cats:
                    ds_categories[CatName(cat_name)] = CatConf(path_in_cfg=("cats",), name=cat_name)
            elif isinstance(raw_ds_cats, dict):
                for cat_name, cat_props in raw_ds_cats.items():
                    if isinstance(cat_props, CatConf):
                        ds_categories[CatName(cat_name)] = cat_props
                    else:
                        if "name" not in cat_props:
                            cat_props = {"name": cat_name, **cat_props}
                        ds_categories[CatName(cat_name)] = CatConf(path_in_cfg=("cats",), **cat_props)
            else:
                raise RunError(f"Invalid category type: {type(raw_ds_cats)}")
            categories_per_ds.append(ds_categories)

        # Allow user to overwrite a given category config for this plot
        if pconf.id and pconf.id in config.custom_plot_config:
            for cat_name, user_cat_props in config.custom_plot_config[pconf.id].items():
                for ds_idx in range(len(categories_per_ds)):
                    if cat_name in categories_per_ds[ds_idx].keys():
                        for prop_name, prop_val in user_cat_props.items():
                            setattr(categories_per_ds[ds_idx][CatName(cat_name)], prop_name, prop_val)

        # Filter data to keep only numerals, remove unknown categories and fill missing with NaNs
        filtered_datasets: List[DatasetT] = []
        for ds_idx, raw_ds in enumerate(raw_datasets):
            filtered_ds: DatasetT = {}
            filtered_datasets.append(filtered_ds)
            for sample_name in list(raw_ds.keys()):
                raw_val_by_cat = raw_ds[sample_name]
                filtered_val_by_cat = {}
                for cat_id, _ in categories_per_ds[ds_idx].items():
                    # Remove categories that are not in the categories_per_ds, and fill missing with NaNs
                    val = raw_val_by_cat.get(cat_id, None)
                    if val is not None and not isinstance(val, (float, int)):
                        # Try to parse
                        try:
                            val = int(val)
                        except ValueError:
                            try:
                                val = float(val)
                            except ValueError:
                                val = None
                    if val is None:
                        val = float("nan")
                    elif isinstance(val, float):
                        if math.floor(val) == val:
                            val = int(val)
                    filtered_val_by_cat[cat_id] = val
                # Remove samples with no data
                if all(math.isnan(v) for v in filtered_val_by_cat.values()):
                    continue
                filtered_datasets[ds_idx][sample_name] = filtered_val_by_cat

        return BarPlotInputData(
            anchor=plot_anchor(pconf),
            pconfig=pconf,
            data=filtered_datasets,
            cats=categories_per_ds,
        )

    def to_df(self):
        """
        Save plot data to a parquet file using a tabular representation that's
        optimized for cross-run analysis.
        """
        # Create a list of records for each sample/category pair
        records = []
        for ds_idx, dataset in enumerate(self.data):
            # Get dataset label if available
            dataset_label = None
            if self.pconfig.data_labels and ds_idx < len(self.pconfig.data_labels):
                label = self.pconfig.data_labels[ds_idx]
                if isinstance(label, dict) and "name" in label:
                    dataset_label = label["name"]
                elif isinstance(label, str):
                    dataset_label = label
                else:
                    dataset_label = f"Dataset {ds_idx + 1}"
            else:
                dataset_label = f"Dataset {ds_idx + 1}"

            for sample_name, sample_data in dataset.items():
                for category_name, value in sample_data.items():
                    # Get category configuration if available
                    cat_conf = None
                    if ds_idx < len(self.cats):
                        cat_conf = self.cats[ds_idx].get(category_name)

                    record = {
                        "anchor": self.anchor,
                        "dataset_idx": ds_idx,
                        "dataset_label": dataset_label,
                        "sample_name": str(sample_name),
                        "category_name": str(category_name),
                        "value": value,
                    }

                    # Store category configuration as JSON if available
                    if cat_conf:
                        record["cat_meta"] = cat_conf.model_dump_json()

                    # Add plot configuration metadata
                    if hasattr(self.pconfig, "title") and self.pconfig.title:
                        record["plot_title"] = self.pconfig.title
                    if hasattr(self.pconfig, "xlab") and self.pconfig.xlab:
                        record["x_label"] = self.pconfig.xlab
                    if hasattr(self.pconfig, "ylab") and self.pconfig.ylab:
                        record["y_label"] = self.pconfig.ylab

                    records.append(record)

        return pd.DataFrame(records)

    @classmethod
    def from_df(cls, df: pd.DataFrame, pconfig: Union[Dict, BarPlotConfig], anchor: Anchor) -> "BarPlotInputData":
        """
        Load plot data from a parquet file.
        """
        import json

        # Handle None case
        if df is None or df.empty:
            if isinstance(pconfig, dict):
                pconf = cast(BarPlotConfig, BarPlotConfig.from_pconfig_dict(pconfig))
            else:
                pconf = cast(BarPlotConfig, pconfig)
            return cls(anchor=anchor, pconfig=pconf, data=[], cats=[])

        if isinstance(pconfig, dict):
            pconf = cast(BarPlotConfig, BarPlotConfig.from_pconfig_dict(pconfig))
        else:
            pconf = cast(BarPlotConfig, pconfig)

        # Group by dataset_idx to rebuild data structure
        datasets: List[DatasetT] = []
        cats_per_dataset: List[Dict[CatName, CatConf]] = []

        max_dataset_idx = df["dataset_idx"].max() if not df.empty else 0

        for ds_idx in range(int(max_dataset_idx) + 1):
            ds_group = df[df["dataset_idx"] == ds_idx] if not df.empty else pd.DataFrame()

            # Skip empty datasets
            if ds_group.empty:
                datasets.append({})
                cats_per_dataset.append({})
                continue

            dataset = {}
            cats_dict = {}

            # Process each sample in this dataset
            for sample_name, sample_group in ds_group.groupby("sample_name"):
                sample_data = {}

                # Process each category for this sample
                for _, row in sample_group.iterrows():
                    category_name = row["category_name"]
                    value = row["value"]
                    sample_data[CatName(str(category_name))] = value

                    # Process category metadata if available
                    if "cat_meta" in row and pd.notna(row["cat_meta"]):
                        try:
                            cat_meta = json.loads(row["cat_meta"])
                            # Use CatConf to ensure proper validation
                            if category_name not in cats_dict:
                                cats_dict[category_name] = CatConf(
                                    name=cat_meta.get("name", category_name),
                                    color=cat_meta.get("color"),
                                    path_in_cfg=("cats",),
                                )
                        except (json.JSONDecodeError, TypeError):
                            # Fallback if JSON parsing fails
                            if category_name not in cats_dict:
                                cats_dict[category_name] = CatConf(name=category_name, path_in_cfg=("cats",))

                # Add sample data if not empty
                if sample_data:
                    dataset[SampleName(str(sample_name))] = sample_data

            datasets.append(dataset)
            cats_per_dataset.append(cats_dict)

        return cls(
            anchor=anchor,
            pconfig=pconf,
            data=datasets,
            cats=cats_per_dataset,
        )

    @classmethod
    def merge(
        cls,
        old_data: "BarPlotInputData",
        new_data: "BarPlotInputData",
    ) -> "BarPlotInputData":
        """
        Merge data from a previous run with the current run.

        This allows for comparing bar plot data across multiple runs.
        """
        # Load dataframes from parquet files instead of using in-memory structures
        old_df = None
        if old_data.anchor in report.plot_input_data:
            old_path = report.plot_input_data[old_data.anchor]
            if Path(old_path).exists():
                try:
                    old_df = pd.read_parquet(old_path)
                except Exception as e:
                    logger.debug(f"Failed to load parquet for {old_data.anchor}: {str(e)}")

        # Create dataframe for new data
        new_df = new_data.to_df()

        # If new_df is empty, return early
        if new_df.empty:
            return old_data

        # Add timestamp and run_id to new data
        new_df["timestamp"] = datetime.now().isoformat()

        if hasattr(config, "kwargs") and "run_id" in config.kwargs:
            new_df["run_id"] = config.kwargs["run_id"]

        # If we have both old and new data, merge them
        merged_df = new_df
        if old_df is not None and not old_df.empty:
            # Make sure old data has run_id
            if "run_id" not in old_df.columns and hasattr(config, "kwargs") and "run_id" in config.kwargs:
                old_df["run_id"] = "previous_run"  # Default value for old data without run_id

            # Make sure old data has timestamp
            if "timestamp" not in old_df.columns:
                old_df["timestamp"] = datetime.now().isoformat()

            # Combine the dataframes, keeping all rows
            merged_df = pd.concat([old_df, new_df], ignore_index=True)

            # For duplicates (same sample, category, dataset), keep the latest version
            if "timestamp" in merged_df.columns:
                # Sort by timestamp (newest last)
                merged_df.sort_values("timestamp", inplace=True)

                # Group by the key identifiers and keep the last entry (newest)
                dedupe_columns = ["dataset_idx", "sample_name", "category_name"]
                if "run_id" in merged_df.columns:
                    dedupe_columns.append("run_id")

                merged_df = merged_df.drop_duplicates(subset=dedupe_columns, keep="last")

        # Save the merged dataframe
        file_path = tmp_dir.parquet_dir() / f"{new_data.anchor}.parquet"
        report.plot_input_data[new_data.anchor] = str(file_path)
        merged_df.to_parquet(file_path, compression="gzip")

        return cls.from_df(
            merged_df,
            new_data.pconfig,
            new_data.anchor,
        )


def plot(
    data: Union[InputDatasetT, Sequence[InputDatasetT]],
    cats: Optional[Union[InputCategoriesT, Sequence[InputCategoriesT]]] = None,
    pconfig: Optional[Union[Dict[str, Any], BarPlotConfig]] = None,
) -> Union["BarPlot", str, None]:
    """
    Create a horizontal bar graph. Also save data to intermediate format.

    Expects a 2D dict of sample data. Also, can take info about categories. There are quite a
    few options of how to use this function, see the docs for details.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
                 Can supply a list of dicts and will have buttons to switch
    :param cats: optional list or dict with plot categories
    :param pconfig: optional dict with config key:value pairs
    :return: HTML and JS, ready to be inserted into the page
    """
    # We want to be permissive to user inputs - but normalizing them now to simplify further processing
    inputs = BarPlotInputData.create(data, cats, pconfig)
    inputs = BarPlotInputData.merge_with_previous(inputs)
    if inputs.is_empty():
        return None

    # Parse the data into a chart friendly format
    scale = mqc_colour.mqc_colour_scale("plot_defaults")  # to add colors to the categories if not set
    plot_samples: List[List[SampleName]] = list()
    plot_data: List[List[CatDataDict]] = list()
    for ds_idx, d in enumerate(inputs.data):
        ordered_samples_names: List[SampleName] = [SampleName(s) for s in d.keys()]
        if isinstance(d, OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif inputs.pconfig.sort_samples:
            ordered_samples_names = natsorted([SampleName(s) for s in d.keys()])
        cat_data_dicts: List[CatDataDict] = list()
        sample_d_count: Dict[SampleName, int] = dict()
        for cat_idx, cat_name in enumerate(inputs.cats[ds_idx].keys()):
            cat_data: List[Union[int, float]] = list()
            cat_count = 0
            for s in ordered_samples_names:
                cat_data.append(d[SampleName(s)][cat_name])
                if s not in sample_d_count:
                    sample_d_count[s] = 0
                cat_count += 1
                sample_d_count[s] += 1
            if cat_count > 0:
                if inputs.pconfig.hide_zero_cats is False or not all(x == 0 for x in cat_data if not math.isnan(x)):
                    color: str = inputs.cats[ds_idx][cat_name].color or scale.get_colour(cat_idx, lighten=1)
                    this_dict: CatDataDict = {
                        "name": inputs.cats[ds_idx][cat_name].name,
                        "color": color,
                        "data": cat_data,
                        "data_pct": [],
                    }
                    cat_data_dicts.append(this_dict)

        # Remove empty samples
        for sample_name, cnt in sample_d_count.items():
            if cnt == 0:
                sample_idx = ordered_samples_names.index(sample_name)
                del ordered_samples_names[sample_idx]
                for cat_data_idx, _ in enumerate(cat_data_dicts):
                    del cat_data_dicts[cat_data_idx]["data"][sample_idx]
        if len(cat_data_dicts) > 0:
            plot_samples.append(ordered_samples_names)
            plot_data.append(cat_data_dicts)

    return BarPlot.create(
        cats_lists=plot_data,
        samples_lists=plot_samples,
        pconfig=inputs.pconfig,
        anchor=inputs.anchor,
    )


class Category(BaseModel):
    name: str
    color: str
    data: List[float]
    data_pct: List[float]


class Dataset(BaseDataset):
    cats: List[Category]
    samples: List[str]

    @staticmethod
    def create(
        dataset: BaseDataset,
        cats: Sequence[CatDataDict],
        samples: Sequence[str],
    ) -> "Dataset":
        # Need to reverse samples as the bar plot will show them reversed
        samples = list(reversed(samples))
        fixed_cats: List[Category] = []
        for input_cat in cats:
            if "name" not in input_cat:
                raise ValueError(f"Bar plot {dataset.plot_id}: missing 'name' key in category")

            # Split long category names
            name = "<br>".join(split_long_string(input_cat["name"]))

            # Reformat color to be ready to add alpha in Plotly-JS
            color = spectra.html(input_cat["color"])
            color_str = ",".join([f"{int(float(x) * 256)}" for x in color.rgb])

            # Reverse the data to match the reversed samples
            cat: Category = Category(
                name=name,
                color=color_str,
                data=list(reversed(input_cat["data"])),
                data_pct=list(reversed(input_cat["data_pct"])) if "data_pct" in input_cat else [],
            )

            # Check that the number of samples is the same for all categories
            assert len(samples) == len(cat.data)
            if cat.data_pct:
                assert len(samples) == len(cat.data_pct)

            fixed_cats.append(cat)

        dataset = Dataset(
            **dataset.model_dump(),
            cats=fixed_cats,
            samples=samples,
        )

        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log=False,
        is_pct=False,
        **kwargs,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        fig = go.Figure(layout=layout)

        for cat in self.cats:
            data = cat.data_pct if is_pct else cat.data

            params = copy.deepcopy(self.trace_params)
            assert cat.color is not None
            params["marker"]["color"] = f"rgb({cat.color})"
            fig.add_trace(
                go.Bar(
                    y=self.samples,
                    x=data,
                    name=cat.name,
                    meta=cat.name,
                    **params,
                ),
            )
        return fig

    def save_data_file(self) -> None:
        val_by_cat_by_sample: Dict[str, Dict[str, str]] = defaultdict(dict)
        for cat in self.cats:
            for d_idx, d_val in enumerate(cat.data):
                s_name = self.samples[d_idx]
                val_by_cat_by_sample[s_name][cat.name] = str(d_val)
        report.write_data_file(val_by_cat_by_sample, self.uid)

    def format_dataset_for_ai_prompt(self, pconfig: PConfig, keep_hidden: bool = True) -> str:
        """Format dataset as a markdown table"""
        prompt = ""
        prompt += "|Sample|" + "|".join(cat.name for cat in self.cats) + "|\n"
        prompt += "|---|" + "|".join("---" for _ in self.cats) + "|\n"

        suffix = ""
        if not pconfig.cpswitch_c_active:
            suffix += "%"
            if self.layout["xaxis"]["ticksuffix"] and self.layout["xaxis"]["ticksuffix"] != "%":
                suffix += " " + self.layout["xaxis"]["ticksuffix"]
        else:
            if self.layout["xaxis"]["ticksuffix"]:
                suffix += " " + self.layout["xaxis"]["ticksuffix"]

        for sidx, sample in enumerate(self.samples):
            presudonym = report.anonymize_sample_name(sample)
            prompt += (
                f"|{presudonym}|"
                + "|".join(
                    self.fmt_value_for_llm(
                        (cat.data if not pconfig.cpswitch or pconfig.cpswitch_c_active else cat.data_pct)[sidx]
                    )
                    + suffix
                    for cat in self.cats
                )
                + "|\n"
            )
        return prompt


class BarPlot(Plot[Dataset, BarPlotConfig]):
    datasets: List[Dataset]

    def samples_names(self) -> List[SampleName]:
        names: List[SampleName] = []
        for ds in self.datasets:
            names.extend(SampleName(sample) for sample in ds.samples)
        return names

    @staticmethod
    def create(
        cats_lists: Sequence[Sequence[CatDataDict]],
        samples_lists: Sequence[Sequence[SampleNameT]],
        pconfig: BarPlotConfig,
        anchor: Anchor,
    ) -> "BarPlot":
        """
        :param cats_lists: each dataset is a list of dicts with the keys: {name, color, data},
            where `name` is the category name, `color` is the color of the bar,
            and `data` is a list of values for each sample. Each outer list will
            correspond a separate tab.
        :param samples_lists: list of lists of bar names (that is, sample names). Similarly,
            each outer list will correspond to a separate tab.
        :param pconfig: Plot configuration dictionary
        """
        if len(cats_lists) != len(samples_lists):
            raise ValueError("Number of datasets and samples lists do not match")

        model: Plot[Dataset, BarPlotConfig] = Plot.initialize(
            plot_type=PlotType.BAR,
            pconfig=pconfig,
            anchor=anchor,
            n_samples_per_dataset=[len(x) for x in samples_lists],
            axis_controlled_by_switches=["xaxis"],
            default_tt_label="%{meta}: <b>%{x}</b>",
            defer_render_if_large=False,  # We hide samples on large bar plots, so no need to defer render
            flat_if_very_large=True,  # However, the data is still embedded into the HTML, and we don't want the report size to inflate
        )

        model.datasets = [
            Dataset.create(d, cats=cats, samples=samples)
            for d, cats, samples in zip(model.datasets, cats_lists, samples_lists)
        ]

        # Set the barmode
        barmode = pconfig.stacking  # stacking, but drawing negative values below zero
        if barmode is None:  # For legacy reasons, interpreting non-default None as "group"
            barmode = "group"  # side by side
        if barmode == "normal":  # Legacy
            barmode = "relative"

        max_n_cats = max([len(dataset.cats) for dataset in model.datasets])

        # Set height to also be proportional to the number of cats to fit a legend
        HEIGHT_PER_LEGEND_ITEM = 19
        legend_height = HEIGHT_PER_LEGEND_ITEM * max_n_cats

        max_n_samples = max(len(x) for x in samples_lists) if len(samples_lists) > 0 else 0
        height = determine_barplot_height(
            max_n_samples=max_n_samples,
            # Group mode puts each category in a separate bar, so need to multiply by the number of categories
            max_bars_in_group=max_n_cats if barmode == "group" else 1,
            legend_height=legend_height,
        )

        model.layout.update(
            height=height,
            barmode=barmode,
            bargroupgap=0,
            bargap=0.2,
            yaxis=dict(
                showgrid=False,
                categoryorder="trace",  # keep sample order
                automargin=True,  # to make sure there is enough space for ticks labels
                title=None,
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
            # Re-initiate legend to reset to default legend location on the top right
            legend=go.layout.Legend(
                # We use legend groups with subplots to simulate standard legend interactivity
                # like we had a standard bar graph without subplots. We need to remove the space
                # between the legend groups to make it look like a single legend.
                tracegroupgap=0,
                # Plotly plots the grouped bar graph in a reversed order in respect to
                # the legend, so reversing the legend to match it:
                traceorder="normal" if barmode != "group" else "reversed",
            ),
            hovermode="y unified",
            hoverlabel=dict(
                bgcolor="rgba(255, 255, 255, 0.8)",
                font=dict(color="black"),
            ),
            showlegend=pconfig.use_legend if pconfig.use_legend is not None else True,
        )

        if getattr(config, "barplot_legend_on_bottom", False):
            model.layout.update(
                legend=go.layout.Legend(
                    orientation="h",
                    x=0.5,
                    xanchor="center",
                    y=-0.5,
                    yanchor="top",
                ),
            )

        for dataset in model.datasets:
            if barmode == "group":
                # max category
                xmax_cnt = max(max(cat.data[i] for cat in dataset.cats) for i in range(len(dataset.samples)))
                xmin_cnt = min(min(cat.data[i] for cat in dataset.cats) for i in range(len(dataset.samples)))
            else:
                # max sum of all categories across all samples
                xmax_cnt = max(
                    sum(cat.data[i] if cat.data[i] > 0 else 0 for cat in dataset.cats)
                    for i in range(len(dataset.samples))
                )
                xmin_cnt = min(
                    sum(cat.data[i] if cat.data[i] < 0 else 0 for cat in dataset.cats)
                    for i in range(len(dataset.samples))
                )

            minallowed = 0 if xmin_cnt > 0 else xmin_cnt  # allow bar to start below zero
            maxallowed = dataset.layout["yaxis"]["autorangeoptions"]["maxallowed"]
            if maxallowed is None:
                maxallowed = xmax_cnt

            dataset.layout.update(
                yaxis=dict(
                    title=None,
                    hoverformat=dataset.layout["xaxis"]["hoverformat"],
                    ticksuffix=dataset.layout["xaxis"]["ticksuffix"],
                ),
                xaxis=dict(
                    title=dict(text=dataset.layout["yaxis"]["title"]["text"]),
                    hoverformat=dataset.layout["yaxis"]["hoverformat"],
                    ticksuffix=dataset.layout["yaxis"]["ticksuffix"],
                    autorangeoptions=dict(
                        minallowed=minallowed,
                        maxallowed=maxallowed,
                    ),
                ),
                showlegend=len(dataset.cats) > 1 if pconfig.use_legend is None else pconfig.use_legend,
            )
            dataset.trace_params.update(
                orientation="h",
                marker=dict(line=dict(width=0)),
                textposition="inside",
                insidetextanchor="start",
            )
            if "hovertemplate" in dataset.trace_params:
                # %{text} doesn't work for unified hovermode:
                dataset.trace_params["hovertemplate"] = dataset.trace_params["hovertemplate"].replace("%{text}", "")

            if dataset.layout["xaxis"]["hoverformat"] is None:
                if all(all(isinstance(x, float) or math.isnan(x) for x in cat.data) for cat in dataset.cats):
                    dataset.layout["xaxis"]["hoverformat"] = ",.2f"
                elif all(all(isinstance(x, int) or math.isnan(x) for x in cat.data) for cat in dataset.cats):
                    dataset.layout["xaxis"]["hoverformat"] = ",.0f"

        # Expand data with zeroes if there are fewer values than samples
        for dataset in model.datasets:
            for cat in dataset.cats:
                if len(cat.data) < len(dataset.samples):
                    cat.data.extend([0] * (len(dataset.samples) - len(cat.data)))

        # Calculate and save percentages
        if model.add_pct_tab:
            for _, dataset in enumerate(model.datasets):
                # Count totals for each category
                sums: List[float] = [0 for _ in dataset.cats[0].data]
                for cat in dataset.cats:
                    for sample_idx, val in enumerate(cat.data):
                        if not math.isnan(val):
                            sums[sample_idx] += abs(val)

                # Now, calculate percentages for each category
                for cat in dataset.cats:
                    values = [x for x in cat.data]
                    for sample_idx, val in enumerate(values):
                        sum_for_sample = sums[sample_idx]
                        if sum_for_sample == 0:
                            values[sample_idx] = 0
                        else:
                            values[sample_idx] = float(val + 0.0) / float(sum_for_sample) * 100.0
                    cat.data_pct = values

                if barmode == "group":
                    # calculating the min percentage range as well because it will be negative for negative values
                    dataset.pct_range["xaxis"]["min"] = min(
                        min(cat.data_pct[i] for cat in dataset.cats) for i in range(len(dataset.samples))
                    )
                else:
                    dataset.pct_range["xaxis"]["min"] = min(
                        sum(cat.data_pct[i] if cat.data_pct[i] < 0 else 0 for cat in dataset.cats)
                        for i in range(len(dataset.samples))
                    )

        if model.add_log_tab:
            # Sorting from small to large so the log switch makes sense
            for dataset in model.datasets:
                dataset.cats.sort(key=lambda cat: sum(cat.data))
                # But reversing the legend so the largest bars are still on the top
                model.layout.legend.traceorder = "reversed"

        return BarPlot(**model.__dict__)

    def _plot_ai_header(self) -> str:
        result = super()._plot_ai_header()
        if self.pconfig.ylab:
            result += f"Values: {self.pconfig.ylab}\n"
        return result
