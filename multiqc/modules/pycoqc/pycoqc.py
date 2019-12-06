#!/usr/bin/env python

""" MultiQC module to parse output from pycoQC """

from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict
import yaml
import logging

log = logging.getLogger(__name__)




class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pycoQC', anchor='pycoqc',
        href="https://a-slide.github.io/pycoQC/",
        info="Computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data")
        self.pycoqc_data = {}
        for f in  self.find_log_files('pycoqc'):
            self.pycoqc_data = self.parse_logs(f['f'])

        if len(self.pycoqc_data) == 0:  # Maybe move this?
            raise UserWarning

        self.pycoqc_info = self.pycoqc_data['pycoqc']

        self.pycoqc_all_reads = self.pycoqc_data['All Reads']
        self.all_reads_run_data = self.pycoqc_all_reads['run']
        self.all_reads_basecall = self.pycoqc_all_reads['basecall']
        self.all_reads_median_length = self.median(self.all_reads_basecall['len_percentiles'])
        self.all_reads_median_qual = self.median(self.all_reads_basecall['qual_score_percentiles'])

        self.pycoqc_passed = self.pycoqc_data['Pass Reads']
        self.passed_run_data = self.pycoqc_passed['run']
        self.passed_basecall = self.pycoqc_passed['basecall']
        self.passed_median_length = self.median(self.passed_basecall['len_percentiles'])
        self.passed_median_qual = self.median(self.passed_basecall['qual_score_percentiles'])

        data = {
           "All Reads" : {
               'first_col': self.all_reads_basecall['reads_number'],
               'second_col': self.all_reads_basecall['bases_number'],
               'third_col': self.all_reads_median_length,
               'fourth_col': self.all_reads_median_qual,
               'fifth_col': self.all_reads_basecall['N50'],
               'sixth_col': self.all_reads_run_data['run_duration'],
               'seventh_col': self.all_reads_run_data['active_channels'],
               'eighth_col': self.all_reads_run_data['runid_number'],
               'ninth_col': self.all_reads_run_data['barcodes_number'],
            },
            "Passed Reads" : {
                'first_col': self.passed_basecall['reads_number'],
                'second_col': self.passed_basecall['bases_number'],
                'third_col': self.passed_median_length,
                'fourth_col': self.passed_median_qual,
                'fifth_col': self.passed_basecall['N50'],
                'sixth_col': self.passed_run_data['run_duration'],
                'seventh_col': self.passed_run_data['active_channels'],
                'eighth_col': self.passed_run_data['runid_number'],
                'ninth_col': self.passed_run_data['barcodes_number'],
            },

        }

        headers = OrderedDict()
        headers['first_col'] = {
            'title': 'Reads',
            'description': 'Number of reads',
            'scale': 'False'
        }
        headers['second_col'] = {
            'title': 'Bases',
            'description': 'Number of bases',
            'scale': 'False'
        }
        headers['third_col'] = {
            'title': 'Median Read Length',
            'description': 'Median Read Length',
            'scale': 'False'
        }
        headers['fourth_col'] = {
            'title': 'Median PHRED score',
            'description': 'Median PHRED score',
            'scale': 'False'
        }
        headers['fifth_col'] = {
            'title': 'N50',
            'description': 'N50',
            'scale': 'False'
        }
        headers['sixth_col'] = {
            'title': 'Run Duration',
            'description': 'Run Duration in hours',
            'scale': 'False'
        }
        headers['seventh_col'] = {
            'title': 'Active Channels',
            'description': 'Number of active channels',
            'scale': 'False'
        }
        headers['eighth_col'] = {
            'title': 'Runids',
            'description': 'Number of runids',
            'scale': 'False'
        }
        headers['ninth_col'] = {
            'title': 'Barcodes',
            'description': 'Number of barcodes',
            'scale': 'False'
        }

        self.general_stats_addcols(data, headers)


    def parse_logs(self, f):
        data = yaml.load(f, Loader=yaml.FullLoader)
        return data

    def median(self, values):
        nr = len(values)
        sorted_values = sorted(values)
        return (sum(sorted_values[nr//2-1:nr//2+1])/2.0, sorted_values[nr//2])[nr % 2] if nr else None
