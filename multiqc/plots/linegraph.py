""" MultiQC functions to plot a linegraph """

import inspect
import logging
import re
from typing import List, Dict, Union, Tuple

from multiqc.utils import config, mqc_colour, report
from multiqc.plots.plotly import line

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(data: Union[List[Dict[str, List]], Dict[str, List]], pconfig=None):
    """
    Plot a line graph with X,Y data.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: optional dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML and JS, ready to be inserted into the page
    """
    # Don't just use {} as the default argument as it's mutable. See:
    # http://python-guide-pt-br.readthedocs.io/en/latest/writing/gotchas/
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

    # Given one dataset - turn it into a list
    if not isinstance(data, list):
        data = [data]
    if "data_labels" in pconfig:
        if len(pconfig["data_labels"]) != len(data):
            raise ValueError(
                f"Length of data_labels does not match the number of datasets. "
                f"Please check your module code and ensure that the data_labels "
                f"list is the same length as the data list: {len(pconfig['data_labels'])} != {len(data)}. "
                f"pconfig={pconfig}"
            )
        pconfig["data_labels"] = [dl if isinstance(dl, dict) else {"name": dl} for dl in pconfig["data_labels"]]

    # Validate config if linting
    if config.strict:
        # Get module name
        modname = ""
        callstack = inspect.stack()
        for n in callstack:
            if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                callpath = n[1].split("multiqc/modules/", 1)[-1]
                modname = f">{callpath}< "
                break
        # Look for essential missing pconfig keys
        for k in ["id", "title", "ylab"]:
            if k not in pconfig and any(k not in dl for dl in pconfig.get("data_labels", [])):
                errmsg = f"LINT: {modname}Linegraph pconfig was missing key '{k}'"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
        # Check plot title format
        if not re.match(r"^[^:]*\S: \S[^:]*$", pconfig.get("title", "")):
            errmsg = "LINT: {} Linegraph title did not match format 'Module: Plot Name' (found '{}')".format(
                modname, pconfig.get("title", "")
            )
            logger.error(errmsg)
            report.lint_errors.append(errmsg)

    # Smooth dataset if requested in config
    if pconfig.get("smooth_points", None) is not None:
        for i, data_by_sample in enumerate(data):
            data[i] = smooth_line_data(data_by_sample, pconfig["smooth_points"])

    datasets: List[List[Dict]] = []
    for ds_idx, data_by_sample in enumerate(data):
        this_list_of_series: List[Dict] = []
        for s in sorted(data_by_sample.keys()):
            # Ensure any overwritten conditionals from data_labels (e.g. ymax) are taken in consideration
            series_config = pconfig.copy()
            pairs: List[Tuple[Union[float, int, str], Union[float, int]]] = []
            maxval = 0
            x_are_categories = pconfig.get("categories", False)
            if "data_labels" in pconfig:
                x_are_categories = pconfig["data_labels"][ds_idx].get("categories", x_are_categories)
            # Discard > ymax or just hide?
            # If it never comes back into the plot, discard. If it goes above then comes back, just hide.
            discard_ymax = None
            discard_ymin = None
            xs = data_by_sample[s].keys()
            if not x_are_categories:
                xs = sorted(xs)

            for k in xs:
                if not x_are_categories:
                    if "xmax" in series_config and float(k) > float(series_config["xmax"]):
                        continue
                    if "xmin" in series_config and float(k) < float(series_config["xmin"]):
                        continue
                if data_by_sample[s][k] is not None and "ymax" in series_config:
                    if float(data_by_sample[s][k]) > float(series_config["ymax"]):
                        discard_ymax = True
                    elif discard_ymax is True:
                        discard_ymax = False
                if data_by_sample[s][k] is not None and "ymin" in series_config:
                    if float(data_by_sample[s][k]) > float(series_config["ymin"]):
                        discard_ymin = True
                    elif discard_ymin is True:
                        discard_ymin = False

            # Build the plot data structure
            for k in xs:
                if not x_are_categories and k is not None:
                    if "xmax" in series_config and float(k) > float(series_config["xmax"]):
                        continue
                    if "xmin" in series_config and float(k) < float(series_config["xmin"]):
                        continue
                if data_by_sample[s][k] is not None:
                    if (
                        "ymax" in series_config
                        and float(data_by_sample[s][k]) > float(series_config["ymax"])
                        and discard_ymax is not False
                    ):
                        continue
                    if (
                        "ymin" in series_config
                        and float(data_by_sample[s][k]) < float(series_config["ymin"])
                        and discard_ymin is not False
                    ):
                        continue
                pairs.append((k, data_by_sample[s][k]))
                try:
                    maxval = max(maxval, data_by_sample[s][k])
                except TypeError:
                    pass
            if maxval > 0 or series_config.get("hide_empty") is not True:
                this_series: Dict = {"name": s, "data": pairs}
                try:
                    this_series["color"] = series_config["colors"][s]
                except KeyError:
                    pass
                this_list_of_series.append(this_series)
        datasets.append(this_list_of_series)

    # Add on annotation data series
    try:
        if pconfig.get("extra_series"):
            extra_series = pconfig["extra_series"]
            if isinstance(pconfig["extra_series"], dict):
                extra_series = [[pconfig["extra_series"]]]
            elif isinstance(pconfig["extra_series"], list) and isinstance(pconfig["extra_series"][0], dict):
                extra_series = [pconfig["extra_series"]]
            for i, es in enumerate(extra_series):
                for s in es:
                    datasets[i].append(s)
    except Exception:
        pass

    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    for si, sd in enumerate(datasets):
        for di, data_by_sample in enumerate(sd):
            data_by_sample.setdefault("color", scale.get_colour(di, lighten=1))

    # Make a plot - template custom, or interactive or flat
    mod = get_template_mod()
    if "linegraph" in mod.__dict__ and callable(mod.linegraph):
        try:
            return mod.linegraph(datasets, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return line.plot(datasets, pconfig)


def smooth_line_data(data: Dict[str, Dict], numpoints: int) -> Dict[str, Dict[int, int]]:
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
    for s_name, d in data.items():
        # Check that we need to smooth this data
        if len(d) <= numpoints or len(d) == 0:
            smoothed_data[s_name] = d
            continue

        binsize = (len(d) - 1) / (numpoints - 1)
        first_element_indices = [round(binsize * i) for i in range(numpoints)]
        smoothed_d = {x: y for i, (x, y) in enumerate(d.items()) if i in first_element_indices}
        smoothed_data[s_name] = smoothed_d

    return smoothed_data
