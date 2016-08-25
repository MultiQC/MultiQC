# -*- coding: utf-8 -*-
from functools import partial

from multiqc import plots
from .utils import find_header

COLUMNS = {
    'snps': {'name': 'SNPs', 'key': 'nSNPs'},
    'mnps': {'name': 'MNPs', 'key': 'nMNPs'},
    'insertions': {'name': 'Insertions', 'key': 'nInsertions'},
    'deletions': {'name': 'Deletions', 'key': 'nDeletions'},
    'complex': {'name': 'Complex', 'key': 'nComplex'},
    'symbolic': {'name': 'Symbolic', 'key': 'nSymbolic'},
    'mixed': {'name': 'Mixed', 'key': 'nMixed'},
    'nocalls': {'name': 'No-calls', 'key': 'nNoCalls'},
}

parse = partial(find_header, '#:GATKTable:CountVariants')


def values(novelties):
    """Parse data per sample for compoverlap output."""
    for novelty in novelties:
        if novelty['Novelty'] == 'all':
            data = {column: int(novelty[data['key']]) for column, data
                    in COLUMNS.items()}
    return data


def plot(data):
    """Create bargraph plot of different variant types."""
    keys = {column: {'name': data['name']} for column, data
            in COLUMNS.items()}
    # Config for the plot
    plot_conf = {
        'id': 'gatk_varianteval_variant_plot',
        'title': 'GATK VariantEval variant types',
        'ylab': '# Reads',
        'cpswitch_counts_label': 'Number of Reads'
    }
    return plots.bargraph.plot(data, keys, plot_conf)
