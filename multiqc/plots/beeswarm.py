#!/usr/bin/env python

""" MultiQC functions to plot a beeswarm group """

import logging
import random
from collections import defaultdict

from multiqc.plots import table_object
from multiqc.utils import config, report, util_functions

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"


def plot(data, headers=None, pconfig=None):
    """Helper HTML for a beeswarm plot.
    :param data: A list of data dicts
    :param headers: A list of Dicts / OrderedDicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :return: HTML string
    """
    if headers is None:
        headers = []
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if "id" in pconfig and pconfig["id"] and pconfig["id"] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig["id"]].items():
            pconfig[k] = v

    # Make a datatable object
    dt = table_object.datatable(data, headers, pconfig)

    return make_plot(dt)


def make_plot(dt):
    bs_id = dt.pconfig.get("id", "table_{}".format("".join(random.sample(letters, 4))))

    # Sanitise plot ID and check for duplicates
    bs_id = report.save_htmlid(bs_id)

    categories = []
    s_names = []
    data = []
    dt.raw_vals = defaultdict(lambda: dict())
    for idx, hs in enumerate(dt.headers):
        for k, header in hs.items():
            bcol = "rgb({})".format(header.get("colour", "204,204,204"))

            categories.append(
                {
                    "namespace": header["namespace"],
                    "title": header["title"],
                    "description": header["description"],
                    "max": header["dmax"],
                    "min": header["dmin"],
                    "suffix": header.get("suffix", ""),
                    "decimalPlaces": header.get("decimalPlaces", "2"),
                    "bordercol": bcol,
                }
            )

            # Add the data
            thisdata = []
            these_snames = []
            for s_name, samp in dt.data[idx].items():
                if k in samp:
                    val = samp[k]
                    dt.raw_vals[s_name][k] = val

                    if "modify" in header and callable(header["modify"]):
                        val = header["modify"](val)

                    thisdata.append(val)
                    these_snames.append(s_name)

            data.append(thisdata)
            s_names.append(these_snames)

    if len(s_names) == 0:
        logger.warning("Tried to make beeswarm plot, but had no data")
        return '<p class="text-danger">Error - was not able to plot data.</p>'

    # Plot HTML
    html = """<div class="hc-plot-wrapper"{height}>
        <div id="{bid}" class="hc-plot not_rendered hc-beeswarm-plot"><small>loading..</small></div>
    </div>""".format(
        bid=bs_id,
        height=f' style="height:{dt.pconfig["height"]}px"' if "height" in dt.pconfig else "",
    )

    report.num_hc_plots += 1

    report.plot_data[bs_id] = {"plot_type": "beeswarm", "samples": s_names, "datasets": data, "categories": categories}

    # Save the raw values to a file if requested
    if dt.pconfig.get("save_file") is True:
        fn = dt.pconfig.get("raw_data_fn", "multiqc_{}".format(bs_id))
        util_functions.write_data_file(dt.raw_vals, fn)
        report.saved_raw_data[fn] = dt.raw_vals

    return html
