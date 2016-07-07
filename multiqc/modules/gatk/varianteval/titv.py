# -*- coding: utf-8 -*-
from functools import partial

from .utils import find_header

parse = partial(find_header, '#:GATKTable:TiTvVariantEvaluator')


def values(novelties):
    """Parse data per sample for TiTvVariantEvaluator output."""
    data = {}
    reference = 'unknown'
    for novelty in novelties:
        if novelty['Novelty'] == 'known':
            reference = novelty['CompRod']
            data['known_titv'] = float(novelty['tiTvRatio'])
        elif novelty['Novelty'] == 'novel':
            data['novel_titv'] = float(novelty['tiTvRatio'])
    return data, reference


def general_headers(reference):
    """Prepare headers for the general overview header."""
    headers = {}
    headers['known_titv'] = {
        'title': 'TiTV ratio (known)',
        'description': "TiTV ratio from variants found in '{}'".format(reference),
        'min': 0,
        'scale': 'Blues',
        'shared_key': 'titv_ratio'
    }
    headers['novel_titv'] = {
        'title': 'TiTV ratio (novel)',
        'description': "TiTV ratio from variants NOT found in '{}'".format(reference),
        'min': 0,
        'scale': 'Blues',
        'shared_key': 'titv_ratio'
    }
    return headers
