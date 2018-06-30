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
        self.fastp_all_data = dict() # to save whole JSON
        for f in self.find_log_files('fastp', filehandles=True):
            self.parse_fastp_log(f)

        # Filter to strip out ignored sample names
        self.fastp_data = self.ignore_samples(self.fastp_data)

        if len(self.fastp_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastp_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.fastp_all_data, 'multiqc_fastp')

        # Basic Stats Table
        self.fastp_general_stats_table()

        # Filtering statistics bar plot
        self.add_section (
           name = 'Filtered Reads',
           anchor = 'fastp-filtered-reads-chart',
           description = 'Filtering statistics of sampled reads.',
           plot = self.fastp_filtered_reads_chart()
        )

    def parse_fastp_log(self, f):
        """ Parse the JSON output from fastp and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except:
            log.warn("Could not parse fastp JSON: '{}'".format(f['fn']))
            return None

        self.add_data_source(f, f['s_name'])
        self.fastp_data[f['s_name']] = {}
        self.fastp_all_data[f['s_name']] = parsed_json

        # Parse filtering_result
        try:
            for k in parsed_json['filtering_result']:
                self.fastp_data[f['s_name']][k] = float(parsed_json['filtering_result'][k])
        except KeyError:
            log.debug("fastp JSON did not have 'filtering_result' key: '{}'".format(f['fn']))

        # Parse duplication
        try:
            self.fastp_data[f['s_name']]['pct_duplication'] = float(parsed_json['duplication']['rate'] * 100.0)
        except KeyError:
            log.debug("fastp JSON did not have a 'duplication' key: '{}'".format(f['fn']))

        # Parse after_filtering
        try:
            for k in parsed_json['summary']['after_filtering']:
                self.fastp_data[f['s_name']][k] = float(parsed_json['summary']['after_filtering'][k])
        except KeyError:
            log.debug("fastp JSON did not have a 'summary'-'after_filtering' keys: '{}'".format(f['fn']))


        # Parse data required to calculate Pct reads surviving
        try:
            self.fastp_data[f['s_name']]['pre_reads'] = float(parsed_json['summary']['before_filtering']['total_reads'])
        except KeyError:
            log.debug("Could not find pre-filtering # reads: '{}'".format(f['fn']))

        try:
            self.fastp_data[f['s_name']]['pct_surviving'] = (float(self.fastp_data[f['s_name']]['total_reads']) / float(self.fastp_data[f['s_name']]['pre_reads'])) * 100.0
        except KeyError:
            log.debug("Could not calculate 'pct_surviving': {}".format(f['fn']))

        # Parse adapter_cutting
        for k in parsed_json['adapter_cutting']:
            try:
                self.fastp_data[f['s_name']][k] = float(parsed_json['adapter_cutting'][k])
            except (ValueError, TypeError):
                pass
            except KeyError:
                log.debug("fastp JSON did not have a 'adapter_cutting' key, skipping: '{}'".format(f['fn']))

        try:
            self.fastp_data[f['s_name']]['pct_adapter'] = (float(self.fastp_data[f['s_name']]['adapter_trimmed_reads']) / float(self.fastp_data[f['s_name']]['pre_reads'])) * 100.0
        except KeyError:
            log.debug("Could not calculate 'pct_adapter': {}".format(f['fn']))


    def fastp_general_stats_table(self):
        """ Take the parsed stats from the fastp report and add it to the
        General Statistics table at the top of the report """

        headers = OrderedDict()
        headers['pct_surviving'] = {
            'title': '% PF',
            'description': 'Percent reads passing filter',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'BuGn',
        }
        headers['gc_content'] = {
            'title': 'GC content',
            'description': 'GC content after filtering',
            'min': 0,
            'scale': 'Blues'
        }
        headers['pct_duplication'] = {
            'title': '% Duplication',
            'description': 'Duplication rate in filtered reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['q30_rate'] = {
            'title': '{} Q30 reads'.format(config.read_count_prefix),
            'description': 'Reads > Q30 after filtering ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count'
        }
        headers['q30_bases'] = {
            'title': '{} Q30 bases'.format(config.base_count_prefix),
            'description': 'Bases > Q30 after filtering ({})'.format(config.base_count_desc),
            'min': 0,
            'modify': lambda x: x * config.base_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'base_count'
        }
        headers['low_quality_reads'] = {
            'title': '{} Low quality'.format(config.read_count_prefix),
            'description': 'Low quality reads ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'RdYlGn-rev',
            'shared_key': 'read_count'
        }
        headers['too_many_N_reads'] = {
            'title': 'Too many N',
            'description': 'Reads with too many N',
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'RdYlGn-rev',
            'shared_key': 'read_count'
        }
        headers['too_short_reads'] = {
            'title': 'Too short',
            'description': 'Reads too short',
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'RdYlGn-rev',
            'shared_key': 'read_count'
        }
        headers['pct_adapter'] = {
            'title': '% Adapter',
            'description': 'Percentage adapter-trimmed reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev',
        }

        self.general_stats_addcols(self.fastp_data, headers)

    def fastp_filtered_reads_chart(self):
       """ Function to generate the fastp filtered reads bar plot """
       # Specify the order of the different possible categories
       keys = OrderedDict()
       keys['passed_filter_reads'] =  { 'name': 'Passed Filter' }
       keys['low_quality_reads'] =    { 'name': 'Low Quality' }
       keys['too_many_N_reads'] =     { 'name': 'Too Many N' }
       keys['too_short_reads'] =      { 'name': 'Too short' }

       # Config for the plot
       pconfig = {
           'id': 'fastp_filtered_reads_plot',
           'title': 'Fastp: Filtered Reads',
           'ylab': '# Reads',
           'cpswitch_counts_label': 'Number of Reads',
           'hide_zero_cats': False,
       }
       return bargraph.plot(self.fastp_data, keys, pconfig)
