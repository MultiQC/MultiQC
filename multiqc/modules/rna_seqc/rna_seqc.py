#!/usr/bin/env python

""" MultiQC module to parse output from RNA-SeQC """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='RNA-SeQC', anchor='rna_seqc',
        href='http://archive.broadinstitute.org/cancer/cga/rna_seqc',
        info="is a java program which computes a series of quality control metrics for RNA-seq data.")

        # Parse metrics information. JSON win!
        self.rna_seqc_metrics = dict()
        for f in self.find_log_files(config.sp['rna_seqc']['metrics']):
            self.parse_metrics(f)

        if len(self.rna_seqc_metrics) == 0 and len(self.rna_seqc_fld) == 0:
            log.debug("Could not find any RNA-SeQC data in {}".format(config.analysis_dir))
            raise UserWarning

        if len(self.rna_seqc_metrics) > 0:
            log.info("Found {} metrics reports".format(len(self.rna_seqc_metrics)))
            self.write_data_file(self.rna_seqc_metrics, 'multiqc_rna_seqc')

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['Exonic Rate'] = {
            'title': '% Exonic',
            'description': 'Exonic rate',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0,
            'format': '{:.1f}%'
        }
        headers['Intronic Rate'] = {
            'title': '% Intronic',
            'description': 'Intronic rate',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0,
            'format': '{:.1f}%'
        }
        headers['Genes Detected'] = {
            'title': '# Genes',
            'description': 'Number of genes detected',
            'min': 0,
            'scale': 'Bu',
            'format': '{:.0f}'
        }
        self.general_stats_addcols(self.rna_seqc_metrics, headers)

    def parse_metrics(self, f):
        """
        Parse the metrics.tsv file from RNA-SeQC
        """
        headers = None
        for l in f['f'].splitlines():
            s = l.split("\t")
            if headers is None:
                headers = s
            else:
                s_name = s[ headers.index('Sample') ]
                data = dict()
                for idx, h in enumerate(headers):
                    data[h] = s[idx]
                self.rna_seqc_metrics[s_name] = data

