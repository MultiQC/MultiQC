#!/usr/bin/env python

""" MultiQC module to parse output from Afterqc """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import json
from distutils.version import StrictVersion

from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Afterqc module class, parses stderr logs.
    Also understands logs saved by Trim Galore!
    (which contain afterqc logs)
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='AfterQC', anchor='afterqc',
        href='http://opengene.org/AfterQC/',
        info="Automatic Filtering, Trimming, Error Removing and Quality Control for fastq data.")

        # Find and load any Afterqc reports
        self.afterqc_data = dict()
        self.afterqc_length_counts = dict()
        self.afterqc_length_exp = dict()
        self.afterqc_length_obsexp = dict()

        for f in self.find_log_files('afterqc', filehandles=True):
            self.parse_afterqc_logs(f)

        # Filter to strip out ignored sample names
        self.afterqc_data = self.ignore_samples(self.afterqc_data)

        if len(self.afterqc_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.afterqc_data)))

        # Write parsed report data to a file
        self.write_data_file(self.afterqc_data, 'multiqc_afterqc')

        # Basic Stats Table
        self.afterqc_general_stats_table()
        # Alignment bar plot
        self.add_section (
            name = 'Bad Reads',
            anchor = 'after_qc',
            plot = self.after_qc_bad_reads_chart()
        )

    def parse_afterqc_logs(self, f):
        """ Go through log file looking for afterqc output """
        fh = json.load(f['f'])

        s_name = f['s_name']

        self.add_data_source(f, s_name)
        self.afterqc_data[s_name] = {}
        self.afterqc_data[s_name]['bad_reads'] = fh['summary']['bad_reads']
        self.afterqc_data[s_name]['bad_reads_with_bad_barcode'] = fh['summary']['bad_reads_with_bad_barcode']
        self.afterqc_data[s_name]['bad_reads_with_bad_overlap'] = fh['summary']['bad_reads_with_bad_overlap']
        self.afterqc_data[s_name]['bad_reads_with_bad_read_length'] = fh['summary']['bad_reads_with_bad_read_length']
        self.afterqc_data[s_name]['bad_reads_with_low_quality'] = fh['summary']['bad_reads_with_low_quality']
        self.afterqc_data[s_name]['bad_reads_with_polyX'] = fh['summary']['bad_reads_with_polyX']
        self.afterqc_data[s_name]['bad_reads_with_reads_in_bubble'] = fh['summary']['bad_reads_with_reads_in_bubble']
        self.afterqc_data[s_name]['bad_reads_with_too_many_N'] = fh['summary']['bad_reads_with_too_many_N']
        self.afterqc_data[s_name]['good_bases'] = fh['summary']['good_bases']
        self.afterqc_data[s_name]['good_reads'] = fh['summary']['good_reads']
        self.afterqc_data[s_name]['readlen'] = fh['summary']['readlen']
        self.afterqc_data[s_name]['total_bases'] = fh['summary']['total_bases']
        self.afterqc_data[s_name]['total_reads'] = fh['summary']['total_reads']


    def afterqc_general_stats_table(self):
        """ Take the parsed stats from the Afterqc report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers["good_bases"] = {
            'title': '{} Good Bases'.format(config.read_count_prefix),
            'description': 'Good Bases ({})'.format(config.read_count_desc),
            'min':0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Pubu'
        }
        headers["total_bases"] = {
            'title': '{} Total Bases'.format(config.read_count_prefix),
            'description': 'Total Bases ({})'.format(config.read_count_desc),
            'min':0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'BuGn'
        }
        headers["good_reads"] = {
            'title': '{} Good Reads'.format(config.read_count_prefix),
            'description': 'Good Reads ({})'.format(config.read_count_desc),
            'min':0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'RdYlGn',

        }
        headers["total_reads"] = {
            'title': '{} Total Reads'.format(config.read_count_prefix),
            'description': 'Total Reads ({})'.format(config.read_count_desc),
            'min':0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues',
            'shared_key': 'read_count'
        }
        headers["readlen"] = {
            'title': 'Read Length',
            'description': 'Read Length',
            'min':0,
            'suffix': ' bp',
            'format': '{:,.0f}',
            'scale': 'YlGn'
        }
        self.general_stats_addcols(self.afterqc_data, headers)

    def after_qc_bad_reads_chart(self):
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['bad_reads'] =      { 'color': '#437bb1', 'name': 'Bad Reads' }
        keys['bad_reads_with_bad_barcode'] =          { 'color': '#7cb5ec', 'name': 'Bad Reads With Bad Barcode' }
        keys["bad_reads_with_bad_overlap"] = { 'color': '#f7a35c', 'name': 'Bad Reads With Bad Overlap'}
        keys["bad_reads_with_bad_read_length"] = { 'color': '#e63491', 'name': 'Bad Reads With Bad Read Length'}
        keys["bad_reads_with_low_quality"] = { 'color': '#b1084c', 'name': 'Bad Reads With Low Quality'}
        keys["bad_reads_with_polyX"] = { 'color': '', 'name': 'Bad Reads With Polyx'}
        keys["bad_reads_with_reads_in_bubble"] = { 'color': '#7f0000', 'name': 'Bad Reads With Reads In Bubble'}
        # Config for the plot
        pconfig = {
            'id': 'after_qc_bad_reads_plot',
            'title': 'After QC Bad Reads',
            'ylab': '# Bad Reads',
            'cpswitch_counts_label': 'Number of Reads',
            'hide_zero_cats': False,
        }
        return bargraph.plot(self.afterqc_data, keys, pconfig)
