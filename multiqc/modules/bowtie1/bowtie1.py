#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 1 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

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
        for f in self.find_log_files(config.sp['bowtie']):
            if f['fn'] in fn_ignore:
                log.debug('Skipping file because looks like tophat log: {}/{}'.format(f['root'], f['fn']))
                continue
            parsed_data = self.parse_bowtie_logs(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.bowtie_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.bowtie_data[f['s_name']] = parsed_data

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
        # Only one section, so add to the intro
        self.intro += self.bowtie_alignment_plot()


    def parse_bowtie_logs(self, s):
        # Check that this isn't actually Bismark using bowtie
        if s.find('bisulfite', 0) >= 0: return None
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
        for k, r in regexes.items():
            match = re.search(r, s)
            if match:
                parsed_data[k] = float(match.group(1).replace(',', ''))
        if len(parsed_data) == 0:
            return None
        return parsed_data


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
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        headers['reads_aligned'] = {
            'title': 'M Aligned',
            'description': 'reads with at least one reported alignment (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
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
            'title': 'Bowtie 1 Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return plots.bargraph.plot(self.bowtie_data, keys, config)
