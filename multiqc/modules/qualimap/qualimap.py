#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
from collections import OrderedDict
import io
import logging
import os

from collections import defaultdict

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Qualimap is really a collection of separate programs:
    BamQC, RNASeq and Counts.. This module is split into separate
    files to reflect this and help with code organisation. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='QualiMap', anchor='qualimap',
        href="http://qualimap.bioinfo.cipf.es/", 
        info="is a platform-independent application to facilitate the quality"\
        " control of alignment sequencing data and its derivatives like"\
        " feature counts.")
        
        # Global dict used by all submodules
        self.general_stats = defaultdict(dict)
        
        # Initialise the BamQC submodule and parse logs
        from . import QM_BamQC
        QM_BamQC.parse_reports(self)
        
        # Initialise the RNASeq submodule and parse logs
        from . import QM_RNASeq
        QM_RNASeq.parse_reports(self)


        if len(self.general_stats) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.general_stats)))

        # Add the submodule results to the report output
        self.sections = list()
        QM_BamQC.report_sections(self)
        QM_RNASeq.report_sections(self)

        # BamQC General stats
        headers = OrderedDict()
        headers['avg_gc'] = {
            'title': 'Avg. GC',
            'description': 'Average GC content',
            'max': 80,
            'min': 20,
            'format': '{:.0f}%'
        }
        headers['median_insert_size'] = {
            'title': 'Insert Size',
            'description': 'Median Insert Size',
            'min': 0,
            'scale': 'PuOr',
            'format': '{:.0f}'
        }
        headers['thirty_x_pc'] = {
            'title': '&ge; 30X',
            'description': 'Fraction of genome with at least 30X coverage',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        headers['median_coverage'] = {
            'title': 'Coverage',
            'description': 'Median coverage',
            'min': 0,
            'scale': 'RdBu'
        }
        headers['mapped_reads'] = {
            'title': 'Aligned',
            'description': 'Number of mapped reads (millions)',
            'min': 0,
            'scale': 'PuBu',
            'shared_key': 'read_count',
            'modify': lambda x: x / 1000000,
        }
        headers['total_reads'] = {
            'title': 'Total Reads',
            'description': 'Number of reads (millions)',
            'min': 0,
            'scale': 'PuBu',
            'shared_key': 'read_count',
            'modify': lambda x: x / 1000000,
        }
        
        
        # RNASeqQC General stats
        headers['5_3_bias'] = {
            'title': "5'-3' bias"
        }
        headers['reads_aligned_genes'] = {
            'title': 'Reads in Genes',
            'description': 'Reads Aligned - Genes (millions)',
            'min': 0,
            'scale': 'PuBu',
            'shared_key': 'read_count',
            'modify': lambda x: x / 1000000,
        }
        headers['reads_aligned'] = {
            'title': 'Aligned',
            'description': 'Reads Aligned (millions)',
            'min': 0,
            'scale': 'RdBu',
            'shared_key': 'read_count',
            'modify': lambda x: x / 1000000,
        }
        
        self.general_stats_addcols(self.general_stats, headers)

    
