#!/usr/bin/env python

""" MultiQC module to parse output from Fastp """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    fastp module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='fastp', anchor='fastp',
        href='https://github.com/OpenGene/fastp',
        info="An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting...)")

        # Find and load any fastp reports
        self.fastp_data = dict()
        for f in self.find_log_files('fastp', filehandles=True):
            self.parse_fastp_log(f)

        # Filter to strip out ignored sample names
        self.fastp_data = self.ignore_samples(self.fastp_data)

        if len(self.fastp_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastp_data)))

        # Write parsed report data to a file
        self.write_data_file(self.fastp_data, 'multiqc_fastp')

        # Basic Stats Table
        self.fastp_general_stats_table()

        # Alignment bar plot
        #self.add_section (
        #    name = 'Bad Reads',
        #    anchor = 'fastp',
        #    description = 'Filtering statistics of sampled reads.',
        #    plot = self.fastp_bad_reads_chart()
        #)

    def parse_fastp_log(self, f):
        """ Parse the JSON output from fastp and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except:
            log.warn("Could not parse fastp JSON: '{}'".format(f['fn']))
            return None

        if 'filtering_result' in parsed_json:
            summaryk = 'filtering_result'
        else:
            log.warn("fastp JSON did not have a 'filtering_result' key, skipping: '{}'".format(f['fn']))
            return None

        s_name = f['s_name']
        self.add_data_source(f, s_name)
        self.fastp_data[s_name] = {}
        for k in parsed_json[summaryk]:
                try:
                    self.fastp_data[s_name][k] = float(parsed_json[summaryk][k])
                except ValueError:
                    self.fastp_data[s_name][k] = parsed_json[summaryk][k]
        #try:
        #    self.fastp_data[s_name]['pct_good_bases'] = (self.fastp_data[s_name]['good_bases'] / self.fastp_data[s_name]['total_bases']) * 100.0
        #except keyerror:
        #    pass

    def fastp_general_stats_table(self):
        """ Take the parsed stats from the fastp report and add it to the
        General Statistics table at the top of the report """

        headers = OrderedDict()
        #headers['pct_good_bases'] = {
        #    'title': '% Good Bases',
        #    'description': 'Percent Good Bases',
        #    'max': 100,
        #    'min': 0,
        #    'suffix': '%',
        #    'scale': 'BuGn',
        #}
        headers['passed_filter_reads'] = {
            'title': '{} PF reads'.format(config.read_count_prefix),
            'description': 'Reads passing filter ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count'
        }
        headers['low_quality_reads'] = {
            'title': '{} Low quality'.format(config.read_count_prefix),
            'description': 'Low quality reads ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues',
            'shared_key': 'read_count'
        }
        headers['too_many_N_reads'] = {
            'title': 'Too many N',
            'description': 'Reads with too many N',
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues',
            'shared_key': 'read_count'
        }
        headers['too_short_reads'] = {
            'title': 'Too short',
            'description': 'Reads too short',
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues',
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.fastp_data, headers)

    #def fastp_qc_bad_reads_chart(self):
    #    """ Function to generate the fastp bad reads bar plot """
    #    # Specify the order of the different possible categories
    #    keys = OrderedDict()
    #    keys['good_reads'] =                     { 'name': 'Good Reads' }
    #    keys['bad_reads_with_bad_barcode'] =     { 'name': 'Bad Barcode' }
    #    keys['bad_reads_with_bad_overlap'] =     { 'name': 'Bad Overlap' }
    #    keys['bad_reads_with_bad_read_length'] = { 'name': 'Bad Read Length' }
    #    keys['bad_reads_with_low_quality'] =     { 'name': 'Low Quality' }
    #    keys['bad_reads_with_polyX'] =           { 'name': 'PolyX' }
    #    keys['bad_reads_with_reads_in_bubble'] = { 'name': 'Reads In Bubble' }
    #    keys['bad_reads_with_too_many_N'] =      { 'name': 'Too many N' }

    #    # Config for the plot
    #    pconfig = {
    #        'id': 'fastp_bad_reads_plot',
    #        'title': 'fastp: Filtered Reads',
    #        'ylab': '# Reads',
    #        'cpswitch_counts_label': 'Number of Reads',
    #        'hide_zero_cats': False,
    #    }
    #    return bargraph.plot(self.fastp_data, keys, pconfig)
