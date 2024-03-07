""" MultiQC functions to plot a box plot """

from typing import List, Dict, Union, OrderedDict

import inspect
import logging
import re

from multiqc.plots.plotly.box import BoxT
from multiqc.utils import config, report
from multiqc.plots.plotly import box

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(list_of_data_by_sample: Union[Dict[str, BoxT], List[Dict[str, BoxT]]], pconfig=None):
    """
    Plot a box plot. Expects either:
    - a dict mapping sample names to data point lists or dicts,
    - a dict mapping sample names to a dict of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    """
    if pconfig is None:
        pconfig = {}

    # Given one dataset - turn it into a list
    if not isinstance(list_of_data_by_sample, list):
        list_of_data_by_sample = [list_of_data_by_sample]

    for i in range(len(list_of_data_by_sample)):
        if isinstance(list_of_data_by_sample[0], OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconfig.get("sort_samples", True):
            samples = sorted(list(list_of_data_by_sample[0].keys()))
            list_of_data_by_sample[i] = {s: list_of_data_by_sample[i][s] for s in samples}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

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
        for k in ["id", "title"]:
            if k not in pconfig:
                errmsg = f"LINT: {modname}Box plot pconfig was missing key '{k}'"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
        # Check plot title format
        if not re.match(r"^[^:]*\S: \S[^:]*$", pconfig.get("title", "")):
            errmsg = "LINT: {} Box plot title did not match format 'Module: Plot Name' (found '{}')".format(
                modname, pconfig.get("title", "")
            )
            logger.error(errmsg)
            report.lint_errors.append(errmsg)

    # Make a plot - custom, interactive or flat
    mod = get_template_mod()
    if "box" in mod.__dict__ and callable(mod.bargraph):
        try:
            return mod.box(list_of_data_by_sample, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return box.plot(list_of_data_by_sample, pconfig)
