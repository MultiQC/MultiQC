#!/usr/bin/env python

""" MultiQC module to parse ATAC-seq pipeline stats """

from __future__ import print_function
from collections import OrderedDict
import logging
import csv

from multiqc import config
from multiqc.plots import bargraph, linegraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    atacseq module class
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='ATAC-seq Stats', anchor='atacseq',
                                            href='https://github.com/epigen/open_pipelines/blob/master/pipelines/atacseq.md',
                                            info="The ATAC-seq pipeline processes ATAC-seq and DNAse-seq data.")
        log.info('Initialized atacseq module')
        # Parse stats for each sample
        self.atacseq_data = dict()
        for f in self.find_log_files(sp_key='atacseq'):
            self.atacseq_data[f['s_name']] = self.parse_atacseq_stats(f['f'])
        log.info('Found stats file for {} ATAC-seq samples'.format(len(self.atacseq_data)))

        # Raise the not found warning
        if len(self.atacseq_data) == 0:
            raise UserWarning

        # Remove ignored samples if there is any
        self.atacseq_data = self.ignore_samples(self.atacseq_data)
        # Add stats to general table
        self.add_atacseq_to_general_stats()

    def parse_atacseq_stats(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split('\t')
            data[s[0]] = s[1]
        return data

    def add_atacseq_to_general_stats(self):
        data = {}
        for sample_name in self.atacseq_data:
            data[sample_name] = {}
            if 'NSC' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['NSC'])
                except:
                    value = None
                data[sample_name]['NSC'] = value
            if 'RSC' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['RSC'])
                except:
                    value = None
                data[sample_name]['RSC'] = value
            if 'peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['peaks'])
                except:
                    value = None
                data[sample_name]['peaks'] = value
            if 'filtered_peaks' in self.atacseq_data[sample_name]:
                try:
                    value = int(self.atacseq_data[sample_name]['filtered_peaks'])
                except:
                    value = None
                data[sample_name]['filtered_peaks'] = value
            if 'frip' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['frip'])
                except:
                    value = None
                data[sample_name]['frip'] = value
            if 'oracle_frip' in self.atacseq_data[sample_name]:
                try:
                    value = float(self.atacseq_data[sample_name]['oracle_frip'])
                except:
                    value = None
                data[sample_name]['oracle_frip'] = value
        headers = OrderedDict()
        headers['peaks'] = {
            'description': 'Number of detected peaks',
            'title': 'Peaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        }
        headers['filtered_peaks'] = {
            'description': 'Number of peaks remaining after filtering',
            'title': 'Filtered\nPeaks',
            'scale': 'Greens',
            'format': '{:,.0f}'
        }
        headers['NSC'] = {
            'description': 'Normalized Strand Cross-correlation Coefficient',
            'title': 'NSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        }
        headers['RSC'] = {
            'description': 'Relative Strand Cross-correlation Coefficient',
            'title': 'RSC',
            'scale': 'Reds',
            'min': 0.0,
            'max': 2.0,
            'format': '{:.2f}'
        }
        headers['frip'] = {
            'description': 'Fraction of Reads in Peaks',
            'title': 'FRiP',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        }
        headers['oracle_frip'] = {
            'description': 'Fraction of Reads in Oracle Peaks',
            'title': 'Oracle FRiP',
            'scale': 'Reds-rev',
            'min': 0.0,
            'max': 1.0,
            'format': '{:.2f}'
        }
        self.general_stats_addcols(data, headers)
