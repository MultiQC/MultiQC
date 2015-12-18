#!/usr/bin/env python

""" MultiQC module to parse output from methylQA """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='methylQA', anchor='methylqa', 
        target='methylQA', href="http://methylqa.sourceforge.net/", 
        info=" - a methylation sequencing data quality assessment tool.")

        # Find and load any methylQA reports
        self.methylqa_data = dict()
        for f in self.find_log_files(contents_match='# reads processed:'):
            parsed_data = self.parse_methylqa_logs(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.methylqa_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.methylqa_data[f['s_name']] = parsed_data

        if len(self.methylqa_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.methylqa_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.methylqa_data, 'multiqc_methylqa1.txt')

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.methylqa_general_stats_table()

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.methylqa_alignment_plot()


    def parse_methylqa_logs(self, s):
        # Check that this isn't actually Bismark using methylqa
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


    def methylqa_general_stats_table(self):
        """ Take the parsed stats from the methylQA report and add it to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['reads_aligned_percentage'] = {
            'title': '% Aligned',
            'description': '% reads with at least one reported alignment',
            'max': 100,
            'min': 0,
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
        self.general_stats_addcols(self.methylqa_data, headers)

    def methylqa_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': 'Aligned' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': 'Multimapped' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }
        
        # Config for the plot
        config = {
            'title': 'methylQA Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.methylqa_data, keys, config)
