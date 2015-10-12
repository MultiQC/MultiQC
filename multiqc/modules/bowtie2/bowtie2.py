#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 2 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import (config, BaseMultiqcModule)

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name="Bowtie 2", anchor="bowtie2", 
        href='http://bowtie-bio.sourceforge.net/bowtie2/', 
        info="is ultrafast and memory-efficient tool for aligning sequencing"\
                " reads to long reference sequences.")

        # Find and load any Bowtie 2 reports
        self.bowtie2_data = dict()
        for f in self.find_log_files(contents_match='reads; of these:'):
            parsed_data = self.parse_bowtie2_logs(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.bowtie2_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.bowtie2_data[f['s_name']] = parsed_data

        if len(self.bowtie2_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bowtie2_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.bowtie2_data, 'multiqc_bowtie2.txt')
        
        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie2_general_stats_table()

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.bowtie2_alignment_plot()


    def parse_bowtie2_logs(self, s):
        # Check that this isn't actually Bismark using bowtie
        if s.find('Using bowtie 2 for aligning with bismark.', 0) >= 0: return None
        parsed_data = {}
        regexes = {
            'reads_processed': r"(\d+) reads; of these:",
            'reads_aligned': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?exactly 1 time",
            'reads_aligned_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?exactly 1 time",
            'not_aligned': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?0 times",
            'not_aligned_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?0 times",
            'multimapped': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?>1 times",
            'multimapped_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?>1 times",
            'overall_aligned_rate': r"([\d\.]+)% overall alignment rate",
        }

        for k, r in regexes.items():
            match = re.search(r, s)
            if match:
                parsed_data[k] = float(match.group(1).replace(',', ''))
            
        if len(parsed_data) == 0: return None
        parsed_data['reads_other'] = parsed_data['reads_processed'] - parsed_data.get('reads_aligned', 0) - parsed_data.get('not_aligned', 0) - parsed_data.get('multimapped', 0)
        return parsed_data


    def bowtie2_general_stats_table(self):
        """ Take the parsed stats from the Bowtie 2 report and add it to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['overall_aligned_rate'] = {
            'title': '% Aligned',
            'description': 'overall alignment rate',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        headers['reads_aligned'] = {
            'title': 'M Aligned',
            'description': 'reads aligned (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.bowtie2_data, headers)

    def bowtie2_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': '1 Alignment' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': '>1 Alignments' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }
        keys['reads_other'] =   { 'color': '#fd0000', 'name': 'Other' }
        
        # Config for the plot
        config = {
            'title': 'Bowtie 2 Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.bowtie2_data, keys, config)
