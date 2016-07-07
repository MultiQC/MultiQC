# -*- coding: utf-8 -*-
from functools import partial

from .utils import find_header

parse = partial(find_header, '#:GATKTable:CountVariants')


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
