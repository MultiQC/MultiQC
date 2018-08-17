#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output from phantompeakqualtools """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='phantompeakqualtools', anchor='phantompeakqualtools',
        href='https://www.encodeproject.org/software/phantompeakqualtools',
        info="computes informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data.")

        # Parse logs
        self.phantompeakqualtools_data = dict()
        for f in self.find_log_files('phantompeakqualtools/out', filehandles=False):
            self.parse_phantompeakqualtools(f)

        # Parse csv for strand correlation plot
        self.correlation_data = dict()
        for f in self.find_log_files('phantompeakqualtools/csv'):
            self.parse_correlation(f)

        # Filter to strip out ignored sample names
        self.phantompeakqualtools_data  =  self.ignore_samples(self.phantompeakqualtools_data)
        self.correlation_data           =  self.ignore_samples(self.correlation_data)

        # Warning when no files are found
        if max(len(self.phantompeakqualtools_data), len(self.correlation_data)) == 0:
            raise UserWarning

        # Log
        log.info("Found {} logs".format(len(self.phantompeakqualtools_data)))
        log.info("Found {} csv file for strand correlation".format(len(self.correlation_data)))

        # Write parsed data to a file
        self.write_data_file(self.phantompeakqualtools_data, 'multiqc_phantompeakqualtools')
        self.write_data_file(self.correlation_data, 'multiqc_correlationdata')

        # Report section
        self.phantompeakqualtools_general_stats()

        if len(self.correlation_data) > 0:
            self.add_section (
                name = 'Strand Shift Correlation Plot',
                anchor = 'strand_shift_correlation',
                plot = self.strand_shift_correlation_plot()
            )

    # Parse spp.out file from phantompeakqualtools
    def parse_phantompeakqualtools(self, f):
        s_name = self.clean_s_name(f['s_name'], f['root'])
        parsed_data = {}
        lines = f['f'].splitlines()
        for l in lines:
            s = l.split("\t")
            parsed_data['Estimated_Fragment_Length_bp'] = int(s[2].split(",")[0])
            parsed_data['NSC'] = float(s[8])
            parsed_data['RSC'] = float(s[9])
        if len(parsed_data) > 0:
            if s_name in self.phantompeakqualtools_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.phantompeakqualtools_data[s_name] = parsed_data

    # Report fragment length, NSC and RSC in general stat table
    def phantompeakqualtools_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['Estimated_Fragment_Length_bp'] = {
            'title': 'Frag Length',
            'description': 'Estimated fragment length (bp)',
            'min': 0,
            'format': '{:,.0f}'
        }
        headers['NSC'] = {
            'title': 'NSC',
            'description': 'Normalized strand cross-correlation',
            'max': 10,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlGn-rev'
        }
        headers['RSC'] = {
            'title': 'RSC',
            'description': 'Relative strand cross-correlation',
            'max': 10,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.phantompeakqualtools_data, headers)

    # Parse spp.csv file from phantompeakqualtools
    def parse_correlation(self, f):
        s_name = self.clean_s_name(f['s_name'], f['root'])
        parsed_data = {}
        reader = csv.DictReader(f['f'])
        for row in reader:
            print(row['x'])
            parsed_data['strand−shift'] = int(row['x'])
            parsed_data['cross-correlation'] = float(row['y'])
            parsed_data['peak_category'] = row['peak']
        if len(parsed_data) > 0:
            if s_name in self.correlation_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.correlation_data[s_name] = parsed_data

    # Strand shift correlation plot
    def strand_shift_correlation_plot(self):
        """ Generate the strand shift correlation plot"""

        data = dict()
        for s_name in self.correlation_data:
            try:
                data[s_name] = {self.correlation_data[s_name]['strand−shift'] : self.correlation_data[s_name]['cross-correlation']}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for strand shift correlation plot')
            return None

        config = {
            'id': 'strand_shift_correlation_plot',
            'title': 'phantompeakqualtools: Strand Shift Correlation Plot',
            'ylab': 'Cross-correlation',
            'xlab': 'Strand shift (bp)',
            'ymin': 0,
            'xmin': 1,
            'xDecimals': False,
            'tt_label': '<b>Strand shift (bp) {point.x}</b>: {point.y} Cross-correlation',
        }

        return linegraph.plot(data, config)
