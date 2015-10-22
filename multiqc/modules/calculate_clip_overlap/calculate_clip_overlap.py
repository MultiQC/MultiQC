#!/usr/bin/env python

""" MultiQC module to parse output from Calculate_Clip_Overlap """

from __future__ import print_function
from collections import OrderedDict
import io
import json
import logging
import os

from multiqc import (config, BaseMultiqcModule)

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name="Calculate_Clip_Overlap", anchor="calculate_clip_overlap", 
        href='http://github.com/avilella/calculate_clip_overlap/', 
        info="calculates the overlap clipping of a given bam file under "\
         "a variety of different Illumina sequencing conditions.")

        # Find and load any Calculate_Clip_Overlap reports
        self.calculate_clip_overlap_data = dict()

        for f in self.find_log_files('clip_overlap.tsv'):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            parsed_data = self.parse_calculate_clip_overlap_logs(f['f'])
#            if parsed_data is not None:
#                if f['s_name'] in self.calculate_clip_overlap_data:
#                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
#                self.calculate_clip_overlap_data[f['s_name']] = parsed_data

        if len(self.calculate_clip_overlap_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.calculate_clip_overlap_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.calculate_clip_overlap_data, 'multiqc_calculate_clip_overlap.txt')
        
        self.sections = list()

        # # Basic Stats Table
        # # Report table is immutable, so just updating it works
        # self.calculate_clip_overlap_general_stats_table()

        # # Alignment Rate Plot
        # # Only one section, so add to the intro
        # self.intro += self.calculate_clip_overlap_alignment_plot()

        
    def parse_calculate_clip_overlap_logs(self, s):
        # Check that this isn't actually Bismark using calculate_clip_overlap
        parsed_data = {}
        regexes = {
            'obs':   r"obs\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
            '2x75':  r"2x75\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
            '2x100': r"2x100\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
            '2x125': r"2x125\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
            '2x150': r"2x150\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
            '2x250': r"2x250\t\d+\t[\d\.]+\d+\t\d+\t\d+\t([\d\.]+)",
        }

        for k, r in regexes.items():
            match = re.search(r, s)
            if match:
                parsed_data[k] = float(match.group(1).replace(',', ''))
            
#        if len(parsed_data) == 0: return None
#        parsed_data['reads_other'] = parsed_data['reads_processed'] - parsed_data.get('reads_aligned', 0) - parsed_data.get('not_aligned', 0) - parsed_data.get('multimapped', 0)
        return parsed_data

    def calculate_clip_overlap_general_stats_table(self):
        """ Take the parsed stats from the Calculate_Clip_Overlap report and add it to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['nonclipped_efficiency_rate'] = {
            'title': '% NonClipped',
            'description': 'overall rate',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        self.general_stats_addcols(self.calculate_clip_overlap_data, headers)

    def calculate_clip_overlap_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': '1 Alignment' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': '>1 Alignments' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }
        keys['reads_other'] =   { 'color': '#fd0000', 'name': 'Other' }
        
        # Config for the plot
        config = {
            'title': 'Calculate_Clip_Overlap Stats',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.calculate_clip_overlap_data, keys, config)
