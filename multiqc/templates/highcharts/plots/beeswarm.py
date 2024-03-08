#!/usr/bin/env python

""" MultiQC functions to plot a beeswarm group """

import logging
import random
from collections import defaultdict

from multiqc.plots import table_object
from multiqc.utils import config, report, util_functions

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


def plot(dt: table_object.DataTable):
    """
    Make a plot - template custom, or interactive or flat
    """
    categories = []
    s_names = []
    data = []
    dt.raw_vals = defaultdict(lambda: dict())
    for idx, hs in enumerate(dt.headers):
        for k, header in hs.items():
            bcol = f"rgb({header.get('colour', '204,204,204')})"

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

    bs_id = dt.pconfig.get("id", f"table_{''.join(random.sample(letters, 4))}")

    # Sanitise plot ID and check for duplicates
    bs_id = report.save_htmlid(bs_id)

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
        fn = dt.pconfig.get("raw_data_fn", f"multiqc_{bs_id}")
        util_functions.write_data_file(dt.raw_vals, fn)
        report.saved_raw_data[fn] = dt.raw_vals

    return html
