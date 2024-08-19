"""MultiQC functions to plot a linegraph"""

import logging
from typing import List, Dict, Union, Tuple, Sequence, Any, Iterable, TypeVar, Generator

from multiqc import config
from multiqc.plots.plotly import line
from multiqc.plots.plotly.line import LinePlotConfig, Series, ValueT, DatasetT, SeriesConf, XToYDictT, KeyT
from multiqc.utils import mqc_colour

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(
    data: Union[DatasetT, Sequence[DatasetT]],
    pconfig: Union[Dict, LinePlotConfig, None] = None,
) -> Union[line.LinePlot, str]:
    """
    Plot a line graph with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    pconf = LinePlotConfig.from_pconfig_dict(pconfig)

    # Given one dataset - turn it into a list
    raw_dataset_list: List[DatasetT]
    if isinstance(data, Sequence):
        raw_dataset_list = list(data)
    else:
        raw_dataset_list = [data]
    del data

    if pconf.data_labels:
        if len(pconf.data_labels) != len(raw_dataset_list):
            raise ValueError(
                f"Length of data_labels does not match the number of datasets. "
                f"Please check your module code and ensure that the data_labels "
                f"list is the same length as the data list: "
                f"{len(pconf.data_labels)} != {len(raw_dataset_list)}. pconfig={pconf}"
            )
        pconf.data_labels = [dl if isinstance(dl, dict) else {"name": dl} for dl in pconf.data_labels]
    else:
        pconf.data_labels = []

    datasets: List[List[Series]] = []
    for ds_idx, raw_data_by_sample in enumerate(raw_dataset_list):
        list_of_series: List[Series] = []
        for s in sorted(raw_data_by_sample.keys()):
            series: Series = _make_series_dict(pconf, ds_idx, s, raw_data_by_sample[s])
            if pconf.hide_empty and not series.pairs:
                continue
            list_of_series.append(series)
        datasets.append(list_of_series)

    # Add on annotation data series
    # noinspection PyBroadException
    if pconf.extra_series:
        pconf_es: Union[Series, List[Series], List[List[Series]]] = pconf.extra_series
        list_of_list_of_series: List[List[Series]]
        if isinstance(pconf_es, list):
            if isinstance(pconf_es[0], list):
                list_of_list_of_series = pconf_es  # type: ignore
            else:
                list_of_list_of_series = [pconf_es for _ in datasets]  # type: ignore
        else:
            list_of_list_of_series = [[pconf_es] for _ in datasets]

        for i, list_of_raw_series in enumerate(list_of_list_of_series):
            assert isinstance(list_of_raw_series, list)
            for series in list_of_raw_series:
                if i < len(datasets):
                    datasets[i].append(series)

    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    for di, series_by_sample in enumerate(datasets):
        for si, series in enumerate(series_by_sample):
            if not series.color:
                series.color = scale.get_colour(si, lighten=1)

    # Make a plot - template custom, or interactive or flat
    mod = get_template_mod()
    if "linegraph" in mod.__dict__ and callable(mod.linegraph):
        # noinspection PyBroadException
        try:
            return mod.linegraph(datasets, pconf)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return line.plot(datasets, pconf)


def _make_series_dict(
    pconfig: LinePlotConfig,
    ds_idx: int,
    s: str,
    y_by_x: XToYDictT[KeyT, ValueT],
) -> Series:
    pairs: List[Tuple[KeyT, ValueT]] = []

    x_are_categories = pconfig.categories
    ymax = pconfig.ymax
    ymin = pconfig.ymin
    xmax = pconfig.xmax
    xmin = pconfig.xmin
    colors = pconfig.colors
    if pconfig.data_labels:
        dl = pconfig.data_labels[ds_idx]
        if isinstance(dl, dict):
            x_are_categories = dl.get("categories", x_are_categories)
            ymax = dl.get("ymax", ymax)
            ymin = dl.get("ymin", ymin)
            xmax = dl.get("xmax", xmax)
            xmin = dl.get("xmin", xmin)
            colors = dl.get("colors", colors)

    # Discard > ymax or just hide?
    # If it never comes back into the plot, discard. If it goes above then comes back, just hide.
    discard_ymax = None
    discard_ymin = None
    xs = [x for x in y_by_x.keys()]
    if not x_are_categories:
        xs = sorted(xs)

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

    return Series(name=s, pairs=pairs, color=colors.get(s))


def smooth_line_data(data_by_sample: DatasetT, numpoints: int) -> DatasetT:
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
    smoothed_data = dict()
    for s_name, d in data_by_sample.items():
        smoothed_data[s_name] = dict(smooth_array(list(d.items()), numpoints))

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

    result = []
    binsize = (len(items) - 1) / (numpoints - 1)
    first_element_indices = {round(binsize * i) for i in range(numpoints)}
    for i, y in enumerate(items):
        if i in first_element_indices:
            result.append(y)
    return result
