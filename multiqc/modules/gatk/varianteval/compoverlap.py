# -*- coding: utf-8 -*-
from functools import partial
from collections import OrderedDict

from multiqc import plots
from .utils import find_header


parse = partial(find_header, '#:GATKTable:CompOverlap')


def values(novelties):
    """Parse data per sample for compoverlap output."""
    data = {}
    for novelty in novelties:
        if novelty['Novelty'] == 'all':
            data['reference'] = novelty['CompRod']
            data['comp_rate'] = float(novelty['compRate'])
            data['concordant_rate'] = float(novelty['concordantRate'])
            data['eval_variants'] = int(novelty['nEvalVariants'])
            data['novel_sites'] = int(novelty['novelSites'])
        elif novelty['Novelty'] == 'known':
            data['known_sites'] = int(novelty['nEvalVariants'])
    return data


def table(data):
    """Bulid a table from the comp overlaps output."""
    headers = OrderedDict()
    headers['comp_rate'] = {
        'title': 'Compare rate',
        'description': 'Ratio of known variants found in the reference set.',
        'min': 0,
        'max': 1,
        'format': '{:.2f}%',
        'scale': 'Blues',
    }
    headers['concordant_rate'] = {
        'title': 'Concordant rate',
        'description': 'Ratio of variants matching alleles in the reference set.',
        'min': 0,
        'max': 1,
        'format': '{:.2f}%',
        'scale': 'Blues',
    }
    headers['eval_variants'] = {
        'title': 'Evaluated variants',
        'description': 'Number of called variants.',
        'min': 0,
        'format': '{:,}',
    }
    headers['known_sites'] = {
        'title': 'Known sites',
        'description': 'Number of known variants.',
        'min': 0,
        'format': '{:,}',
    }
    headers['novel_sites'] = {
        'title': 'Novel sites',
        'description': 'Number of novel variants.',
        'min': 0,
        'format': '{:,}',
    }
    table_html = plots.table.plot(data, headers)
    return table_html
