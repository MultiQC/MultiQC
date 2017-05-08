#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 1 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Bowtie 1 module, parses stderr logs. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Bowtie 1', anchor='bowtie1',
        target='Bowtie 1', href="http://bowtie-bio.sourceforge.net/",
        info="is an ultrafast, memory-efficient short read aligner.")

        # Find and load any Bowtie reports
        self.bowtie_data = dict()
        fn_ignore = [
            # Tophat log files
            'bowtie.left_kept_reads.log',
            'bowtie.left_kept_reads.m2g_um.log',
            'bowtie.left_kept_reads.m2g_um_seg1.log',
            'bowtie.left_kept_reads.m2g_um_seg2.log',
            'bowtie.right_kept_reads.log',
            'bowtie.right_kept_reads.m2g_um.log',
            'bowtie.right_kept_reads.m2g_um_seg1.log',
            'bowtie.right_kept_reads.m2g_um_seg2.log'
        ]
        for f in self.find_log_files('bowtie'):
            if f['fn'] in fn_ignore:
                log.debug('Skipping file because looks like tophat log: {}/{}'.format(f['root'], f['fn']))
                continue
            # Check that this isn't actually Bismark using bowtie
            if f['f'].find('bisulfite', 0) >= 0:
                log.debug('Skipping file because looks like Bismark log: {}/{}'.format(f['root'], f['fn']))
                continue
            self.parse_bowtie_logs(f)

        # Filter to strip out ignored sample names
        self.bowtie_data = self.ignore_samples(self.bowtie_data)

        if len(self.bowtie_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bowtie_data)))

        # Write parsed report data to a file
        self.write_data_file(self.bowtie_data, 'multiqc_bowtie1')

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie_general_stats_table()

        # Alignment Rate Plot
        self.bowtie_alignment_plot()


    def parse_bowtie_logs(self, f):
        s_name = f['s_name']
        parsed_data = {}
        regexes = {
            'reads_processed': r"# reads processed:\s+(\d+)",
            'reads_aligned': r"# reads with at least one reported alignment:\s+(\d+)",
            'reads_aligned_percentage': r"# reads with at least one reported alignment:\s+\d+\s+\(([\d\.]+)%\)",
            'not_aligned': r"# reads that failed to align:\s+(\d+)",
            'not_aligned_percentage': r"# reads that failed to align:\s+\d+\s+\(([\d\.]+)%\)",
            'multimapped': r"# reads with alignments suppressed due to -m:\s+(\d+)",
            'multimapped_percentage': r"# reads with alignments suppressed due to -m:\s+\d+\s+\(([\d\.]+)%\)"
        }
        for l in f['f'].splitlines():
            # Attempt in vain to find original bowtie1 command, logged by another program
            if 'bowtie' in l and 'q.gz' in l:
                fqmatch = re.search(r"([^\s,]+\.f(ast)?q.gz)", l)
                if fqmatch:
                    s_name = self.clean_s_name(fqmatch.group(1), f['root'])
                    log.debug("Found a bowtie command, updating sample name to '{}'".format(s_name))

            # End of log, reset in case there is another in this file
            if 'Overall time:' in l:
                if len(parsed_data) > 0:
                    if s_name in self.bowtie_data:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.bowtie_data[s_name] = parsed_data
                s_name = f['s_name']
                parsed_data = {}

            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = float(match.group(1).replace(',', ''))

        if len(parsed_data) > 0:
            if s_name in self.bowtie_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.bowtie_data[s_name] = parsed_data


    def bowtie_general_stats_table(self):
        """ Take the parsed stats from the Bowtie report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['reads_aligned_percentage'] = {
            'title': '% Aligned',
            'description': '% reads with at least one reported alignment',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['reads_aligned'] = {
            'title': '{} Aligned'.format(config.read_count_prefix),
            'description': 'reads with at least one reported alignment ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.bowtie_data, headers)

    def bowtie_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': 'Aligned' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': 'Multimapped' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }

        # Config for the plot
        config = {
            'id': 'bowtie1_alignment',
            'title': 'Bowtie 1 Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section( plot = bargraph.plot(self.bowtie_data, keys, config) )
