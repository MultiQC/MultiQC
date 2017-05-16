#!/usr/bin/env python

""" MultiQC functions to plot a beeswarm group """

import logging
import random

from multiqc.utils import report
from multiqc.plots import table_object

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=None, pconfig=None):
    """ Helper HTML for a beeswarm plot.
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

    # Make a datatable object
    dt = table_object.datatable(data, headers, pconfig)

    return make_plot( dt )


def make_plot(dt):

    bs_id = dt.pconfig.get('id', 'table_{}'.format(''.join(random.sample(letters, 4))) )

    # Sanitise plot ID and check for duplicates
    bs_id = report.save_htmlid(bs_id)

    categories = []
    s_names = []
    data = []
    for idx, hs in enumerate(dt.headers):
        for k, header in hs.items():

            bcol = 'rgb({})'.format(header.get('colour', '204,204,204'))

            categories.append({
                'namespace': header['namespace'],
                'title': header['title'],
                'description': header['description'],
                'max': header['dmax'],
                'min': header['dmin'],
                'suffix': header.get('suffix', ''),
                'decimalPlaces': header.get('decimalPlaces', '2'),
                'bordercol': bcol
            });

            # Add the data
            thisdata = []
            these_snames = []
            for (s_name, samp) in dt.data[idx].items():
                if k in samp:

                    val = samp[k]

                    if 'modify' in header and callable(header['modify']):
                        val = header['modify'](val)

                    thisdata.append(val)
                    these_snames.append(s_name)

            data.append(thisdata)
            s_names.append(these_snames)

    # Plot HTML
    html = """<div class="hc-plot-wrapper">
        <div id="{bid}" class="hc-plot not_rendered hc-beeswarm-plot"><small>loading..</small></div>
    </div>""".format(bid=bs_id)

    report.num_hc_plots += 1

    report.plot_data[bs_id] = {
        'plot_type': 'beeswarm',
        'samples': s_names,
        'datasets': data,
        'categories': categories
    }

    return html


