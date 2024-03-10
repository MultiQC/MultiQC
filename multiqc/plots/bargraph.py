""" MultiQC functions to plot a bargraph """


import inspect
import logging
from collections import OrderedDict

import math
import re

from multiqc.utils import config, mqc_colour, report
from multiqc.plots.plotly import bar

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


def plot(data, cats=None, pconfig=None):
    """Plot a horizontal bar graph. Expects a 2D dict of sample
    data. Also, can take info about categories. There are quite a
    few variants of how to use this function, see the docs for details.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
                 Can supply a list of dicts and will have buttons to switch
    :param cats: optional list or dict with plot categories
    :param pconfig: optional dict with config key:value pairs
    :return: HTML and JS, ready to be inserted into the page
    """

    if pconfig is None:
        pconfig = {}

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
        for k in ["id", "title", "ylab"]:
            if k not in pconfig:
                errmsg = f"LINT: {modname}Bargraph pconfig was missing key '{k}'"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
        # Check plot title format
        if not re.match(r"^[^:]*\S: \S[^:]*$", pconfig.get("title", "")):
            errmsg = "LINT: {} Bargraph title did not match format 'Module: Plot Name' (found '{}')".format(
                modname, pconfig.get("title", "")
            )
            logger.error(errmsg)
            report.lint_errors.append(errmsg)

    # Given one dataset - turn it into a list
    if not isinstance(data, list):
        data = [data]

    # Make list of cats from different inputs
    if cats is None:
        cats = list()
    elif not isinstance(cats, list):
        cats = [cats]
    elif cats and isinstance(cats[0], str):
        cats = [cats]
    # Generate default categories if not supplied
    for idx in range(len(data)):
        try:
            cats[idx]
        except IndexError:
            cats.append(list())
            for s in data[idx].keys():
                for k in data[idx][s].keys():
                    if k not in cats[idx]:
                        cats[idx].append(k)

    # If we have cats in lists, turn them into dicts
    for idx, cat in enumerate(cats):
        if isinstance(cat, list):
            newcats = dict()
            for c in cat:
                newcats[c] = {"name": c}
            cats[idx] = newcats
        else:
            for c in cat:
                if "name" not in cat[c]:
                    cats[idx][c]["name"] = c

    # Allow user to overwrite a given category config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            for idx in range(len(cats)):
                if k in cats[idx].keys():
                    for kk, vv in v.items():
                        cats[idx][k][kk] = vv

    # Parse the data into a chart friendly format
    plotsamples = list()
    plotdata = list()
    for idx, d in enumerate(data):
        hc_samples = list(d.keys())
        if isinstance(d, OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconfig.get("sort_samples", True):
            hc_samples = sorted(list(d.keys()))
        hc_data = list()
        sample_dcount = dict()
        for c in cats[idx].keys():
            thisdata = list()
            catcount = 0
            for s in hc_samples:
                if s not in sample_dcount:
                    sample_dcount[s] = 0

                if s not in d or c not in d[s]:
                    # Pad with NaNs when we have missing categories in a sample
                    thisdata.append(float("nan"))
                    continue
                val = d[s][c]
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
                    thisdata.append(float("nan"))
                    continue
                if isinstance(val, float):
                    if math.floor(val) == val:
                        val = int(val)
                thisdata.append(val)
                catcount += 1
                sample_dcount[s] += 1
            if catcount > 0:
                if pconfig.get("hide_zero_cats", True) is False or max(x for x in thisdata if not math.isnan(x)) > 0:
                    thisdict = {"name": cats[idx][c]["name"], "data": thisdata}
                    if "color" in cats[idx][c]:
                        thisdict["color"] = cats[idx][c]["color"]
                    hc_data.append(thisdict)

        # Remove empty samples
        for s, c in sample_dcount.items():
            if c == 0:
                idx = hc_samples.index(s)
                del hc_samples[idx]
                for j, d in enumerate(hc_data):
                    del hc_data[j]["data"][idx]
        if len(hc_data) > 0:
            plotsamples.append(hc_samples)
            plotdata.append(hc_data)

    if len(plotdata) == 0:
        logger.warning(f"Tried to make bar plot, but had no data: {pconfig.get('id')}")
        return '<p class="text-danger">Error - was not able to plot data.</p>'

    # Add colors to the categories if not set. Since the "plot_defaults" scale is
    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    for si, sd in enumerate(plotdata):
        for di, d in enumerate(sd):
            d.setdefault("color", scale.get_colour(di, lighten=1))

    # Make a plot - custom, interactive or flat
    mod = get_template_mod()
    if "bargraph" in mod.__dict__ and callable(mod.bargraph):
        try:
            return mod.bargraph(plotdata, plotsamples, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return bar.plot(plotdata, plotsamples, pconfig)
