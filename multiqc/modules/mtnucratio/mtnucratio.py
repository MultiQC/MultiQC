#!/usr/bin/env python

""" MultiQC module to parse output from mtnucratio """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
import json

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ mtnucratio module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='DeDup', anchor='dedup',
        href="http://www.github.com/apeltzer/DeDup",
        info="is a tool for duplicate removal for merged/collapsed reads in ancient DNA analysis.")

        # Find and load any DeDup reports
        self.dedup_data = dict()

        for f in self.find_log_files('dedup',filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.dedup_data = self.ignore_samples(self.dedup_data)

        if len(self.dedup_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.dedup_data)))

        # Write parsed report data to a file
        self.write_data_file(self.dedup_data, 'multiqc_dedup')

        # Basic Stats Table
        self.dedup_general_stats_table()

        # Alignment Rate Plot
        self.add_section(
            description = 'This plot shows read categories that were either not removed (unique reads) or removed (duplicates).',
            plot = self.dedup_alignment_plot()
        )

    #Parse our nice little JSON file
    def parseJSON(self, f):
        """ Parse the JSON output from mtnucratio and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except Exception as e:
            print(e)
            log.warn("Could not parse mtnucratio JSON: '{}'".format(f['fn']))
            return None

        #Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json['metadata']['sample_name'],'')
        self.add_data_source(f, s_name)

        metrics_dict = parsed_json['metrics']

        #Add all in the main data_table
        self.dedup_data[s_name] = metrics_dict


    def dedup_general_stats_table(self):
        """ Take the parsed stats from the mtnucratio report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['mt_cov_avg'] = {
            'title': 'Average mitochondrial genome coverage',
            'description': 'Average coverage (X) on mitochondrial genome.',
            'min': 0,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.2f}',
        }
        headers['nuc_cov_avg'] = {
            'title': 'Average nuclear genome coverage',
            'description': 'Average coverage (X) on nuclear genome.',
            'min': 1,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.2f}',
        }
        headers['mt_nuc_ratio'] = {
            'title': 'MTNUC Ratio',
            'description': 'Mitochondrial to nuclear reads ratio (MTNUC)',
            'min': 1,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.2f}',
        }
        headers['nucreads'] = {
            'title': 'Reads on nuclear genome',
            'description': 'Reads on the nuclear genome.',
            'min': 1,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.0f}',
        }
        headers['mtreads'] = {
            'title': 'Reads on mitochondrial genome',
            'description': 'Reads on the mitochondrial genome.',
            'min': 1,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.2f}',
        }

        self.general_stats_addcols(self.dedup_data, headers)

    