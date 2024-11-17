"""MultiQC functions to plot a bargraph"""

import logging
import math
from collections import OrderedDict
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union, cast

from importlib_metadata import EntryPoint

from multiqc import config
from multiqc.core.exceptions import RunError
from multiqc.plots.plotly import bar
from multiqc.plots.plotly.bar import BarPlotConfig, CatDataDict
from multiqc.types import SampleName
from multiqc.utils import mqc_colour
from multiqc.validation import ValidatedConfig

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod: Optional[EntryPoint] = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


SampleNameT = Union[SampleName, str]
InputDatasetT = Mapping[SampleNameT, Mapping[str, Union[int, float]]]
DatasetT = Dict[SampleNameT, Dict[str, Union[int, float]]]


class CatConf(ValidatedConfig):
    name: str
    color: Optional[str] = None


# Either a list of strings, or a dictionary mapping category names to their properties dicts or objects
CategoriesT = Union[Sequence[str], Mapping[str, Union[Mapping[str, str], CatConf]]]


def plot(
    data: Union[InputDatasetT, Sequence[InputDatasetT]],
    cats: Optional[Union[CategoriesT, Sequence[CategoriesT]]] = None,
    pconfig: Optional[Union[Dict[str, Any], BarPlotConfig]] = None,
) -> Union[bar.BarPlot, str]:
    """Plot a horizontal bar graph. Expects a 2D dict of sample
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
                ds_categories[cat_name] = CatConf(name=cat_name)
        elif isinstance(raw_ds_cats, dict):
            for cat_name, cat_props in raw_ds_cats.items():
                if isinstance(cat_props, CatConf):
                    ds_categories[cat_name] = cat_props
                else:
                    if "name" not in cat_props:
                        cat_props["name"] = cat_name
                    ds_categories[cat_name] = CatConf(**cat_props, _clss=[bar.BarPlot])
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

    # Parse the data into a chart friendly format
    # To add colors to the categories if not set:
    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    plot_samples: List[List[SampleName]] = list()
    plot_data: List[List[CatDataDict]] = list()
    for ds_idx, d in enumerate(raw_datasets):
        hc_samples: List[SampleName] = [SampleName(s) for s in d.keys()]
        if isinstance(d, OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconf.sort_samples:
            hc_samples = sorted([SampleName(s) for s in d.keys()])
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
                if pconf.hide_zero_cats is False or max(x for x in this_data if not math.isnan(x)) > 0:
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

    # Make a plot - custom, interactive or flat
    mod = get_template_mod()
    if "bargraph" in mod.__dict__ and callable(mod.__dict__["bargraph"]):
        try:
            return mod.__dict__["bargraph"](plot_data, plot_samples, pconfig)
        except:  # noqa: E722s
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return bar.plot(plot_data, plot_samples, pconf)
