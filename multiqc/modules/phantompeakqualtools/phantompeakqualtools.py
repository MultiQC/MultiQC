#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output from phantompeakqualtools """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

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
        for f in self.find_log_files('phantompeakqualtools', filehandles=True):
            self.parse_phantompeakqualtools(f)

        # Filter to strip out ignored sample names
        self.phantompeakqualtools_data = self.ignore_samples(self.phantompeakqualtools_data)

        if len(self.phantompeakqualtools_data) == 0:
            raise UserWarning

        log.info("Found {} logs".format(len(self.phantompeakqualtools_data)))
        self.write_data_file(self.phantompeakqualtools_data, 'multiqc_phantompeakqualtools')

        self.phantompeakqualtools_general_stats()


    def parse_phantompeakqualtools(self, f):
        s_name = f['s_name']
        parsed_data = {}
        lines = f['f'].splitlines()
        for line in lines:
            s = l.split("\t")
            parse_data['Estimated Fragment Length (bp)'] = int(s[2].split(",")[0])
            parse_data['NSC'] = float(s[8])
            parse_data['RSC'] = float(s[9])
        if len(parsed_data) > 0:
            if s_name in self.phantompeakqualtools_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.phantompeakqualtools_data[s_name] = parsed_data


    def phantompeakqualtools_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['est_length'] = {
            'title': 'Estimated Fragment Length (bp)',
            'min': 0,
            'format': '{:,.0f}'
        }
        headers['nsc'] = {
            'title': 'NSC,
            'description': 'Normalized strand cross-correlation',
            'max': 1,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlBu-rev'
        }
        headers['rsc'] = {
            'title': 'RSC',
            'description': 'Relative strand cross-correlation',
            'max': 1,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.phantompeakqualtools_data, headers)
