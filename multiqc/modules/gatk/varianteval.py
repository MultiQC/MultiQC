#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK varianteval """

import logging
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)

class VariantEvalMixin():

    def parse_gatk_varianteval(self):
        """ Find GATK varianteval logs and parse their data """

        self.gatk_varianteval = dict()
        for f in self.find_log_files('gatk/varianteval', filehandles=True):
            parsed_data = parse_single_report(f['f'])
            if len(parsed_data) > 0:
                if f['s_name'] in self.gatk_varianteval:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='varianteval')
                self.gatk_varianteval[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.gatk_varianteval = self.ignore_samples(self.gatk_varianteval)

        if len(self.gatk_varianteval) > 0:

            # Write parsed report data to a file (restructure first)
            self.write_data_file(self.gatk_varianteval, 'multiqc_gatk_varianteval')

            # Get consensus TiTv references
            titv_ref = None
            for s_name in self.gatk_varianteval:
                if titv_ref is None:
                    titv_ref = self.gatk_varianteval[s_name]['titv_reference']
                elif titv_ref != self.gatk_varianteval[s_name]['titv_reference']:
                    titv_ref = 'Multiple'
                    break

            # General Stats Table
            varianteval_headers = dict()
            varianteval_headers['known_titv'] = {
                'title': 'TiTV ratio (known)',
                'description': "TiTV ratio from variants found in '{}'".format(titv_ref),
                'min': 0,
                'scale': 'Blues',
                'shared_key': 'titv_ratio'
            }
            varianteval_headers['novel_titv'] = {
                'title': 'TiTV ratio (novel)',
                'description': "TiTV ratio from variants NOT found in '{}'".format(titv_ref),
                'min': 0,
                'scale': 'Blues',
                'shared_key': 'titv_ratio'
            }
            self.general_stats_addcols(self.gatk_varianteval, varianteval_headers, 'GATK VariantEval')

            # Variant Counts plot
            self.add_section (
                name = 'Variant Counts',
                anchor = 'gatk-count-variants',
                plot = count_variants_barplot(self.gatk_varianteval)
            )

            # Compare Overlap Table
            self.add_section (
                name = 'Compare Overlap',
                anchor = 'gatk-compare-overlap',
                plot = comp_overlap_table(self.gatk_varianteval)
            )


        # Return the number of logs that were found
        return len(self.gatk_varianteval)


def parse_single_report(f):
    """ Parse a gatk varianteval varianteval """

    data = dict()
    in_CompOverlap = False
    in_CountVariants = False
    in_TiTv = False
    for l in f:
        # Detect section headers
        if '#:GATKTable:CompOverlap' in l:
            in_CompOverlap = True
        elif '#:GATKTable:CountVariants' in l:
            in_CountVariants = True
        elif '#:GATKTable:TiTvVariantEvaluator' in l:
            in_TiTv = True
        else:
            # Parse contents using nested loops
            if in_CompOverlap:
                headers = l.split()
                while in_CompOverlap:
                    l = f.readline().strip("\n")
                    d = dict()
                    try:
                        for i, s in enumerate(l.split()):
                            d[headers[i]] = s
                        if d['Novelty'] == 'all':
                            data['reference'] = d['CompRod']
                            data['comp_rate'] = float(d['compRate'])
                            data['concordant_rate'] = float(d['concordantRate'])
                            data['eval_variants'] = int(d['nEvalVariants'])
                            data['novel_sites'] = int(d['novelSites'])
                        elif d['Novelty'] == 'known':
                            data['known_sites'] = int(d['nEvalVariants'])
                    except KeyError:
                        in_CompOverlap = False
            elif in_CountVariants:
                headers = l.split()
                while in_CountVariants:
                    l = f.readline().strip("\n")
                    d = dict()
                    try:
                        for i, s in enumerate(l.split()):
                            d[headers[i]] = s
                        if d['Novelty'] == 'all':
                            data['snps'] = int(d['nSNPs'])
                            data['mnps'] = int(d['nMNPs'])
                            data['insertions'] = int(d['nInsertions'])
                            data['deletions'] = int(d['nDeletions'])
                            data['complex'] = int(d['nComplex'])
                            data['symbolic'] = int(d['nSymbolic'])
                            data['mixed'] = int(d['nMixed'])
                            data['nocalls'] = int(d['nNoCalls'])
                    except KeyError:
                        in_CountVariants = False
            elif in_TiTv:
                headers = l.split()
                data['titv_reference'] = 'unknown'
                while in_TiTv:
                    l = f.readline().strip("\n")
                    d = dict()
                    try:
                        for i, s in enumerate(l.split()):
                            d[headers[i]] = s
                        if d['Novelty'] == 'known':
                            data['titv_reference'] = d['CompRod']
                            data['known_titv'] = float(d['tiTvRatio'])
                        elif d['Novelty'] == 'novel':
                            data['novel_titv'] = float(d['tiTvRatio'])
                    except KeyError:
                        in_TiTv = False

    return data


def count_variants_barplot(data):
    """ Return HTML for the Variant Counts barplot """
    keys = OrderedDict()
    keys['snps'] = {'name': 'SNPs'}
    keys['mnps'] = {'name': 'MNPs'}
    keys['insertions'] = {'name': 'Insertions'}
    keys['deletions'] = {'name': 'Deletions'}
    keys['complex'] = {'name': 'Complex'}
    keys['symbolic'] = {'name': 'Symbolic'}
    keys['mixed'] = {'name': 'Mixed'}
    keys['nocalls'] = {'name': 'No-calls'}

    plot_conf = {
        'id': 'gatk_varianteval_variant_plot',
        'title': 'GATK VariantEval Variant Counts',
        'ylab': '# Variants',
        'cpswitch_counts_label': 'Number of Variants'
    }
    return bargraph.plot(data, keys, plot_conf)


def comp_overlap_table(data):
    """Build a table from the comp overlaps output."""
    headers = OrderedDict()
    headers['comp_rate'] = {
        'title': 'Compare rate',
        'description': 'Ratio of known variants found in the reference set.',
        'namespace': 'GATK',
        'min': 0,
        'max': 100,
        'suffix': '%',
        'format': '{:,.2f}',
        'scale': 'Blues',
    }
    headers['concordant_rate'] = {
        'title': 'Concordant rate',
        'description': 'Ratio of variants matching alleles in the reference set.',
        'namespace': 'GATK',
        'min': 0,
        'max': 100,
        'suffix': '%',
        'format': '{:,.2f}',
        'scale': 'Blues',
    }
    headers['eval_variants'] = {
        'title': 'M Evaluated variants',
        'description': 'Number of called variants (millions)',
        'namespace': 'GATK',
        'min': 0,
        'modify': lambda x: float(x) / 1000000.0
    }
    headers['known_sites'] = {
        'title': 'M Known sites',
        'description': 'Number of known variants (millions)',
        'namespace': 'GATK',
        'min': 0,
        'modify': lambda x: float(x) / 1000000.0
    }
    headers['novel_sites'] = {
        'title': 'M Novel sites',
        'description': 'Number of novel variants (millions)',
        'namespace': 'GATK',
        'min': 0,
        'modify': lambda x: float(x) / 1000000.0
    }
    table_html = table.plot(data, headers, {'id': 'gatk_compare_overlap', 'table_title': 'GATK - Compare Overlap'})
    return table_html
