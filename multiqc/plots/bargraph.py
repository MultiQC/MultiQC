"""MultiQC functions to plot a bargraph"""

import copy
import logging
import math
from collections import OrderedDict, defaultdict
from typing import Any, Dict, List, Literal, Mapping, Optional, Sequence, Tuple, TypedDict, Union, cast

import plotly.graph_objects as go  # type: ignore
import spectra  # type: ignore
from importlib_metadata import EntryPoint
from natsort import natsorted
from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.core.exceptions import RunError
from multiqc.plots.plotly import determine_barplot_height
from multiqc.plots.plotly.plot import (
    BaseDataset,
    PConfig,
    Plot,
    PlotType,
    split_long_string,
)
from multiqc.types import Anchor, SampleName
from multiqc.utils import mqc_colour
from multiqc.validation import ValidatedConfig

logger = logging.getLogger(__name__)

################################################################################
## Input types

SampleNameT = Union[SampleName, str]
InputDatasetT = Mapping[SampleNameT, Mapping[str, Union[int, float]]]
DatasetT = Dict[SampleNameT, Dict[str, Union[int, float]]]


class CatConf(ValidatedConfig):
    name: str
    color: Optional[str] = None

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("cats",), **data)


# Either a list of strings, or a dictionary mapping category names to their properties dicts or objects
CategoriesT = Union[Sequence[str], Mapping[str, Union[Mapping[str, str], CatConf]]]


################################################################################
# Output types


class BarPlotConfig(PConfig):
    stacking: Union[Literal["group", "overlay", "relative", "normal"], None] = "relative"
    hide_empty: Optional[bool] = Field(None, deprecated="hide_zero_cats")
    hide_zero_cats: bool = True
    sort_samples: bool = True
    use_legend: bool = True
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


SampleNameT = Union[SampleName, str]


class CatDataDict(TypedDict):
    name: str
    color: str
    data: List[float]
    data_pct: List[float]


class NormalizedBarPlotData(BaseModel):
    data: List[Dict[SampleName, Dict[str, Union[int, float]]]]
    cats: List[Dict[str, CatConf]]


def plot(
    data: Union[InputDatasetT, Sequence[InputDatasetT]],
    cats: Optional[Union[CategoriesT, Sequence[CategoriesT]]] = None,
    pconfig: Optional[Union[Dict[str, Any], BarPlotConfig]] = None,
) -> Union["BarPlot", str]:
    """
    Plot a horizontal bar graph. Expects a 2D dict of sample
    data. Also, can take info about categories. There are quite a
    few variants of how to use this function, see the docs for details.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
                 Can supply a list of dicts and will have buttons to switch
    :param cats: optional list or dict with plot categories
    :param pconfig: optional dict with config key:value pairs
    :return: HTML and JS, ready to be inserted into the page
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
    raw_cats_per_ds: List[CategoriesT]
    if not cats:
        # Not supplied, generate default categories
        raw_cats_per_ds = []
        for val_by_cat_by_sample in raw_datasets:
            ds_cats: List[str] = []
            for sample_name, val_by_cat in val_by_cat_by_sample.items():
                for cat_name in val_by_cat.keys():
                    if cat_name not in raw_cats_per_ds:
                        ds_cats.append(cat_name)
            raw_cats_per_ds.append(ds_cats)
    elif isinstance(cats, List) and isinstance(cats[0], str):
        # ["Cat1", "Cat2"] - list of strings for one dataset
        raw_cats_per_ds = [[cat_name for cat_name in cast(List[str], cats)]]
    elif isinstance(cats, Sequence):
        # [["Cat1", "Cat2"], {"Cat3": {}, "Cat4": {}}] - list of lists or dicts for multiple datasets
        raw_cats_per_ds = [ds_cats for ds_cats in cast(List[Dict], cats)]
    else:
        raw_cats_per_ds = [cats]

    if len(raw_datasets) > 1 and len(raw_cats_per_ds) == 1:
        raw_cats_per_ds = raw_cats_per_ds * len(raw_datasets)
    elif len(raw_datasets) != len(raw_cats_per_ds):
        raise RunError(
            f"Bar graph: number of dataset and category lists must match, got {len(raw_datasets)} "
            f"datasets and {len(raw_cats_per_ds)} category lists: {raw_cats_per_ds}"
        )

    # Parse the categories into pydantic objects
    categories_per_ds: List[Dict[str, CatConf]] = []
    for raw_ds_cats in raw_cats_per_ds:
        ds_categories: Dict[str, CatConf] = dict()
        if isinstance(raw_ds_cats, list):
            for cat_name in raw_ds_cats:
                ds_categories[cat_name] = CatConf(path_in_cfg=("cats",), name=cat_name)
        elif isinstance(raw_ds_cats, dict):
            for cat_name, cat_props in raw_ds_cats.items():
                if isinstance(cat_props, CatConf):
                    ds_categories[cat_name] = cat_props
                else:
                    if "name" not in cat_props:
                        cat_props = {"name": cat_name, **cat_props}
                    ds_categories[cat_name] = CatConf(path_in_cfg=("cats",), **cat_props)
        else:
            raise RunError(f"Invalid category type: {type(raw_ds_cats)}")
        categories_per_ds.append(ds_categories)

    # Allow user to overwrite a given category config for this plot
    if pconf.id and pconf.id in config.custom_plot_config:
        for cat_name, user_cat_props in config.custom_plot_config[pconf.id].items():
            for ds_idx in range(len(categories_per_ds)):
                if cat_name in categories_per_ds[ds_idx].keys():
                    for prop_name, prop_val in user_cat_props.items():
                        setattr(categories_per_ds[ds_idx][cat_name], prop_name, prop_val)

    pid = pconf.anchor or pconf.id
    if prev_plot := report.plot_by_id.get(Anchor(pid)):
        assert isinstance(prev_plot, BarPlot)
        return BarPlot.update(prev_plot, pconf, categories_per_ds, raw_datasets)

    # Parse the data into a chart friendly format
    scale = mqc_colour.mqc_colour_scale("plot_defaults")  # to add colors to the categories if not set
    plot_samples: List[List[SampleName]] = list()
    plot_data: List[List[CatDataDict]] = list()
    for ds_idx, d in enumerate(raw_datasets):
        hc_samples: List[SampleName] = [SampleName(s) for s in d.keys()]
        if isinstance(d, OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconf.sort_samples:
            hc_samples = natsorted([SampleName(s) for s in d.keys()])
        hc_data: List[CatDataDict] = list()
        sample_d_count: Dict[SampleName, int] = dict()
        for cat_idx, c in enumerate(categories_per_ds[ds_idx].keys()):
            this_data: List[Union[int, float]] = list()
            cat_count = 0
            for s in hc_samples:
                if s not in sample_d_count:
                    sample_d_count[s] = 0

                if s not in d or c not in d[s]:
                    # Pad with NaNs when we have missing categories in a sample
                    this_data.append(float("nan"))
                    continue
                val = d[SampleName(s)][c]
                if not isinstance(val, (float, int)):
                    try:
                        val = int(val)
                    except ValueError:
                        try:
                            val = float(val)
                        except ValueError:
                            val = None
                if val is None:
                    # Pad with NaNs when we have missing categories in a sample
                    this_data.append(float("nan"))
                    continue
                if isinstance(val, float):
                    if math.floor(val) == val:
                        val = int(val)
                this_data.append(val)
                cat_count += 1
                sample_d_count[s] += 1
            if cat_count > 0:
                if pconf.hide_zero_cats is False or not all(x == 0 for x in this_data if not math.isnan(x)):
                    color: str = categories_per_ds[ds_idx][c].color or scale.get_colour(cat_idx, lighten=1)
                    this_dict: CatDataDict = {
                        "name": categories_per_ds[ds_idx][c].name,
                        "color": color,
                        "data": this_data,
                        "data_pct": [],
                    }
                    hc_data.append(this_dict)

        # Remove empty samples
        for sample_name, cnt in sample_d_count.items():
            if cnt == 0:
                sample_idx = hc_samples.index(sample_name)
                del hc_samples[sample_idx]
                for hc_data_idx, _ in enumerate(hc_data):
                    del hc_data[hc_data_idx]["data"][sample_idx]
        if len(hc_data) > 0:
            plot_samples.append(hc_samples)
            plot_data.append(hc_data)

    if len(plot_data) == 0:
        logger.warning(f"Tried to make bar plot, but had no data: {pconf.id}")
        return '<p class="text-danger">Error - was not able to plot data.</p>'

    return BarPlot.create(
        pconfig=pconf,
        cats_lists=plot_data,
        samples_lists=plot_samples,
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
    def update(
        dataset: "Dataset",
        new_pconfig: BarPlotConfig,
        new_cats: Sequence[CatDataDict],
        new_samples: Sequence[SampleNameT],
    ) -> "Dataset":
        """
        Update a dataset with new data. Creates a new instance.
        """
        # Convert plot data back to input form
        data_by_sample_by_cat: Dict[SampleName, Dict[str, Union[int, float]]] = defaultdict(dict)
        cats: List[CatConf] = []

        # Convert previous dataset back to input form
        for cat_obj in dataset.cats:
            cats.append(CatConf(name=cat_obj.name, color=cat_obj.color))
            for sample, val in zip(dataset.samples, cat_obj.data):
                data_by_sample_by_cat[SampleName(sample)][cat_obj.name] = val

        for new_cat in new_cats:
            existing_cat = next((existing_cat for existing_cat in cats if existing_cat.name == new_cat["name"]), None)
            if not existing_cat:
                existing_cat = CatConf(name=new_cat["name"], color=new_cat["color"])
                for new_sample, new_val in zip(samples, new_cat["data"]):
                    existing_cat.data[new_sample] = new_val
            else:
                cat_confs.append(CatConf(name=new_cat["name"], color=new_cat["color"]))
                for sample, val in zip(samples, new_cat["data"]):
                    data_by_sample_by_cat[sample][new_cat["name"]] = val

        combined_samples = dataset.samples[:]
        for new_sample in samples:
            if new_sample not in combined_samples:
                combined_samples.append(new_sample)
        if pconfig.sort_samples:
            combined_samples = natsorted(combined_samples)

        sample_to_cat_to_val: Dict[SampleName, Dict[str, float]] = defaultdict(dict)
        sample_to_cat_to_pct: Dict[SampleName, Dict[str, float]] = defaultdict(dict)
        for cat in dataset.cats:
            for sample, val in zip(dataset.samples, cat.data):
                sample_to_cat_to_val[sample][cat.name] = val

        combined_cats = dataset.cats
        for new_cat in cats:
            # Find if this category already exists
            if matching_old_cat := next((cat for cat in dataset.cats if cat.name == new_cat["name"]), None):
                # Update existing category data with new values
                # Create a mapping of sample to value for the new data
                new_data_by_sample = {sample: value for sample, value in zip(samples, new_cat["data"])}
                for sample_idx_in_old, sample in enumerate(dataset.samples):
                    if sample in new_data_by_sample:
                        matching_old_cat.data[sample_idx_in_old] = new_data_by_sample[sample]

                # Update percentages if they exist
                if new_cat.get("data_pct") and matching_old_cat.data_pct:
                    new_pct_by_sample = {sample: value for sample, value in zip(samples, new_cat["data_pct"])}
                    for sample_dx_in_old, sample in enumerate(dataset.samples):
                        if sample in new_pct_by_sample:
                            matching_old_cat.data_pct[sample_dx_in_old] = new_pct_by_sample[sample]
                combined_cats.append(matching_old_cat)
            else:
                # Add new category
                for prev_sample in dataset.samples:
                    if prev_sample in samples:
                        this_sample_idx_in_new = samples.index(prev_sample)
                        combined_data.append(new_cat["data"][this_sample_idx_in_new])
                        if "data_pct" in new_cat:
                            combined_data_pct.append(new_cat["data_pct"][this_sample_idx_in_new])
                    else:
                        combined_data.append(0)
                        if "data_pct" in new_cat:
                            combined_data_pct.append(0)

                combined_cats.append(
                    Category(
                        name=new_cat["name"],
                        color=new_cat["color"],
                        data=combined_data,
                        data_pct=combined_data_pct if "data_pct" in new_cat else [],
                    )
                )

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


class BarPlot(Plot[Dataset, BarPlotConfig]):
    datasets: List[Dataset]

    @staticmethod
    def update(
        plot: "BarPlot",
        data: Sequence[InputDatasetT],
        cats: Sequence[CatConf],
        pconfig: BarPlotConfig,
    ) -> "BarPlot":
        """
        Create plot from existing instance and new data. Returns a new instance.
        """
        datasets: List[Dataset] = []

        # For each dataset
        for ds_idx, (dataset, ds_data, ds_cats) in enumerate(zip(plot.datasets, data, cats)):
            datasets.append(Dataset.update(dataset, pconfig, ds_cats, ds_data))

        return BarPlot.create(
            pconfig,
            [dataset.cats for dataset in datasets],
            [dataset.samples for dataset in datasets],
        )

        #     # Process each new category
        #     for new_cat in new_cats:
        #         # Find if this category already exists
        #         existing_cat = next((cat for cat in dataset.cats if cat.name == new_cat["name"]), None)
        #         if existing_cat:
        #             # Update existing category data
        #             # Create a mapping of sample to value for the new data
        #             new_data_by_sample = {sample: value for sample, value in zip(new_samples, new_cat["data"])}

        #             # Update the existing data array
        #             for sample_idx, sample in enumerate(dataset.samples):
        #                 if sample in new_data_by_sample:
        #                     existing_cat.data[sample_idx] = new_data_by_sample[sample]

        #             # Update percentages if they exist
        #             if new_cat.get("data_pct") and existing_cat.data_pct:
        #                 new_pct_by_sample = {sample: value for sample, value in zip(new_samples, new_cat["data_pct"])}
        #                 for sample_idx, sample in enumerate(dataset.samples):
        #                     if sample in new_pct_by_sample:
        #                         existing_cat.data_pct[sample_idx] = new_pct_by_sample[sample]
        #             combined_cat_lists.append(existing_cat)
        #         else:
        #             # Add new category
        #             # Pad the data array with zeros for existing samples that aren't in the new data
        #             full_data = []
        #             full_data_pct = []

        #             for sample in dataset.samples:
        #                 if sample in new_samples:
        #                     idx = new_samples.index(sample)
        #                     full_data.append(new_cat["data"][idx])
        #                     if "data_pct" in new_cat:
        #                         full_data_pct.append(new_cat["data_pct"][idx])
        #                 else:
        #                     full_data.append(0)
        #                     if "data_pct" in new_cat:
        #                         full_data_pct.append(0)

        #             combined_cat_lists.append(
        #                 Category(
        #                     name=new_cat["name"],
        #                     color=new_cat["color"],
        #                     data=full_data,
        #                     data_pct=full_data_pct if "data_pct" in new_cat else [],
        #                 )
        #             )

        # return BarPlot.create(pconfig, combined_cat_lists, combined_sample_lists)

    @staticmethod
    def create(
        pconfig: BarPlotConfig,
        cats_lists: Sequence[Sequence[CatDataDict]],
        samples_lists: Sequence[Sequence[SampleNameT]],
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
            showlegend=pconfig.use_legend,
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
                showlegend=len(dataset.cats) > 1,
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
