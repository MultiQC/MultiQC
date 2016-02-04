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

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Qualimap RNASeq reports and parse their data """
    
    self.qualimap_rnaseq_gorigin = dict()
    self.qualimap_rnaseq_cov_hist = dict()
    
    sp = config.sp['qualimap']['rnaseq']
    
    # Find QualiMap reports
    for directory in config.analysis_dir:
        for root, dirnames, filenames in os.walk(directory, followlinks=True):
            raw_data_dir = sp['raw_data']
            for d in dirnames:
                if raw_data_dir in d:
                    raw_data_dir = d
            if sp['rnaseq_results'] in filenames and raw_data_dir in dirnames:
                with io.open(os.path.join(root, sp['rnaseq_results']), 'r') as gr:
                    for l in gr:
                        rhs = l.split(' = ')[-1].replace(',','')
                        num = rhs.split('(')[0]
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root)
                            if s_name in self.qualimap_rnaseq_gorigin:
                                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                            self.qualimap_rnaseq_gorigin[s_name] = dict()
                        
                        # Reads alignment section
                        if 'reads aligned' in l and '(left/right)' not in l:
                            self.general_stats[s_name]['reads_aligned'] = float( num )
                        if 'read pairs aligned' in l:
                            self.general_stats[s_name]['reads_aligned'] = float( num )
                        if 'aligned to genes' in l:
                            self.general_stats[s_name]['reads_aligned_genes'] = float( num )
                        if "5'-3' bias" in l:
                            self.general_stats[s_name]['5_3_bias'] = float( num )
                        
                        # Reads genomic origin section
                        if 'exonic' in l:
                            self.qualimap_rnaseq_gorigin[s_name]['reads_aligned_exonic'] = float( num )
                        if 'intronic' in l:
                            self.qualimap_rnaseq_gorigin[s_name]['reads_aligned_intronic'] = float( num )
                        if 'intergenic' in l:
                            self.qualimap_rnaseq_gorigin[s_name]['reads_aligned_intergenic'] = float( num )
                        if 'overlapping exon' in l:
                            self.qualimap_rnaseq_gorigin[s_name]['reads_aligned_overlapping_exon'] = float( num )
                
                #### Coverage profile
                cov_report = os.path.join(root, raw_data_dir, sp['coverage'])
                if os.path.exists(cov_report):
                    self.qualimap_rnaseq_cov_hist[s_name] = {}
                    with io.open(cov_report, 'r') as fh:
                        next(fh)
                        for l in fh:
                            coverage, count = l.split(None, 1)
                            coverage = int(round(float(coverage)))
                            count = float(count)
                            self.qualimap_rnaseq_cov_hist[s_name][coverage] = count
                        
                    


def report_sections(self):
    """ Add results from Qualimap RNASeq parsing to the report """
    
    # Genomic Origin Bar Graph
    if len(self.qualimap_rnaseq_gorigin) > 0:
        gorigin_cats = OrderedDict()
        gorigin_cats['reads_aligned_intergenic'] = {'name': 'Intergenic'}
        gorigin_cats['reads_aligned_intronic'] = {'name': 'Intronic'}
        gorigin_cats['reads_aligned_exonic'] = {'name': 'Exonic'}
        gorigin_pconfig = {
            'title': 'Genomic Origin',
            'cpswitch_c_active': False
        }
        self.sections.append({
            'name': 'Reads genomic origin',
            'anchor': 'qualimap-reads-genomic-origin',
            'content': self.plot_bargraph(self.qualimap_rnaseq_gorigin, gorigin_cats, gorigin_pconfig)
        })
    
    if len(self.qualimap_rnaseq_cov_hist) > 0:
        self.sections.append({
            'name': 'Coverage Profile Along Genes (total)',
            'anchor': 'qualimap-genome-fraction-coverage',
            'content': self.plot_xy_data(self.qualimap_rnaseq_cov_hist, {
                'title': 'Coverage Profile Along Genes (total)',
                'ylab': 'Coverage',
                'xlab': 'Transcript Position (bp)',
                'ymin': 0,
                'xmin': 0,
                'xmax': 100,
                'tt_label': '<b>{point.x} bp</b>: {point.y:.0f}%',
            })
        })

