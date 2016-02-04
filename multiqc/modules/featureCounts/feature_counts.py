#!/usr/bin/env python

""" MultiQC module to parse output from featureCounts """

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
        super(MultiqcModule, self).__init__(name='featureCounts', 
        anchor='featurecounts', target='Subread featureCounts', 
        href='http://bioinf.wehi.edu.au/featureCounts/', 
        info="is a highly efficient general-purpose read summarization program"\
        " that counts mapped reads for genomic features such as genes, exons,"\
        " promoter, gene bodies, genomic bins and chromosomal locations.")

        # Find and load any featureCounts reports
        self.featurecounts_data = dict()
        for f in self.find_log_files(config['sp']['featurecounts']):
            parsed_data = self.parse_featurecounts_report(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.featurecounts_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.featurecounts_data[f['s_name']] = parsed_data

        if len(self.featurecounts_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.featurecounts_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.featurecounts_data, 'multiqc_featureCounts.txt')

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.featurecounts_stats_table()

        # Assignment bar plot
        # Only one section, so add to the intro
        self.intro += self.featureCounts_chart()



    def parse_featurecounts_report (self, raw_data):
        """ Parse the featureCounts log file. """

        regexes = {
            'Assigned': r"Assigned\s+(\d+)",
            'Unassigned_Ambiguity': r"Unassigned_Ambiguity\s+(\d+)",
            'Unassigned_MultiMapping': r"Unassigned_MultiMapping\s+(\d+)",
            'Unassigned_NoFeatures': r"Unassigned_NoFeatures\s+(\d+)",
            'Unassigned_Unmapped': r"Unassigned_Unmapped\s+(\d+)",
            'Unassigned_MappingQuality': r"Unassigned_MappingQuality\s+(\d+)",
            'Unassigned_FragmentLength': r"Unassigned_Frage?mentLength\s+(\d+)",
            'Unassigned_Chimera': r"Unassigned_Chimera\s+(\d+)",
            'Unassigned_Secondary': r"Unassigned_Secondary\s+(\d+)",
            'Unassigned_Nonjunction': r"Unassigned_Nonjunction\s+(\d+)",
            'Unassigned_Duplicate': r"Unassigned_Duplicate\s+(\d+)",
        }
        parsed_data = {}
        total_count = 0
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
                total_count += float(r_search.group(1))
        parsed_data['Total'] = total_count
        if 'Assigned' in parsed_data:
            parsed_data['percent_assigned'] = (parsed_data['Assigned']/parsed_data['Total']) * 100
        if len(parsed_data) == 0: return None
        return parsed_data


    def featurecounts_stats_table(self):
        """ Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['percent_assigned'] = {
            'title': '% Assigned',
            'description': '% Assigned reads',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        headers['Assigned'] = {
            'title': 'M Assigned',
            'description': 'Assigned reads (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.featurecounts_data, headers)


    def featureCounts_chart (self):
        """ Make the featureCounts assignment rates plot """
        
        # Specify the order of the different possible categories
        keys = ['Assigned', 'Unassigned_Ambiguity', 'Unassigned_MultiMapping', 'Unassigned_NoFeatures',
        'Unassigned_Unmapped', 'Unassigned_MappingQuality', 'Unassigned_FragmentLength', 'Unassigned_Chimera',
        'Unassigned_Secondary', 'Unassigned_Nonjunction', 'Unassigned_Duplicate']
        
        # Config for the plot
        config = {
            'id': 'featureCounts_assignment_plot',
            'title': 'featureCounts Assignments',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.featurecounts_data, keys, config)
