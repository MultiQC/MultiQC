#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Submodule to handle code for Qualimap BamQC """

from __future__ import print_function
from collections import OrderedDict
import io
import logging
import os

from collections import defaultdict

from multiqc import config, BaseMultiqcModule

def parse_reports(self):
    """ Find Qualimap RNASeq reports and parse their data """
    
    # Find QualiMap reports
    for directory in config.analysis_dir:
        for root, dirnames, filenames in os.walk(directory, followlinks=True):
            raw_data_dir = 'raw_data'
            for d in dirnames:
                if raw_data_dir in d:
                    raw_data_dir = d
            if 'rnaseq_qc_results.txt' in filenames and raw_data_dir in dirnames:
                with io.open(os.path.join(root, 'rnaseq_qc_results.txt'), 'r') as gr:
                    for l in gr:
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root)
                            if s_name in self.reads_aligned:
                                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                        if 'reads aligned' in l:
                            self.general_stats[s_name]['reads_aligned'] = float( l.split(' = ')[-1].replace(',','') )
                        if 'aligned to genes' in l:
                            self.general_stats[s_name]['reads_aligned_genes'] = float( l.split(' = ')[-1].replace(',','') )
                        if 'exonic' in l:
                            self.general_stats[s_name]['reads_aligned_exonic'] = float( l.split(' = ')[-1].replace(',','') )
                        if 'intronic' in l:
                            self.general_stats[s_name]['reads_aligned_intronic'] = float( l.split(' = ')[-1].replace(',','') )
                        if 'intergenic' in l:
                            self.general_stats[s_name]['reads_aligned_intergenic'] = float( l.split(' = ')[-1].replace(',','') )
                    


def report_sections(self):
    """ Add results from Qualimap RNASeq parsing to the report """
    # Append to self.sections list
    
    return None



def stats_table(self):
    """ Take the parsed stats from the QualiMap RNASeq report and add them to the
    basic stats table at the top of the report """
    
    headers = OrderedDict()
    headers['reads_aligned'] = {
        'title': 'Aligned',
        'description': 'Reads Aligned (millions)',
        'min': 0,
        'scale': 'RdBu',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    headers['reads_aligned_genes'] = {
        'title': 'A: Genes',
        'description': 'Reads Aligned: Genes (millions)',
        'min': 0,
        'scale': 'PuOr',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    headers['reads_aligned_exonic'] = {
        'title': 'A: Exonic',
        'description': 'Reads Aligned: Exonic (millions)',
        'min': 0,
        'scale': 'RdYlGn',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    headers['reads_aligned_intronic'] = {
        'title': 'A: Intronic',
        'description': 'Reads Aligned: Intronic (millions)',
        'min': 0,
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    headers['reads_aligned_intergenic'] = {
        'title': 'A: Intergenic',
        'description': 'Reads Aligned: Intergenic (millions)',
        'min': 0,
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    self.general_stats_addcols(self.general_stats, headers)

