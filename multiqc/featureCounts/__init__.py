#!/usr/bin/env python

""" MultiQC module to parse output from featureCounts """

from __future__ import print_function
import logging
import re

import multiqc
from multiqc import config

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('featureCounts'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(log)

        # Static variables
        self.name = "featureCounts"
        self.anchor = "featurecounts"
        self.intro = '<p><a href="http://bioinf.wehi.edu.au/featureCounts/" target="_blank">Subread featureCounts</a> \
             is a highly efficient general-purpose read summarization program that counts mapped reads for genomic \
             features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations.</p>'

        # Find and load any featureCounts reports
        self.featurecounts_data = dict()
        for f in self.find_log_files('_counts.txt.summary'):
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
            'Unassigned_FragementLength': r"Unassigned_FragementLength\s+(\d+)",
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
        if len(parsed_data) == 0: return None
        return parsed_data


    def featurecounts_stats_table(self):
        """ Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report """

        config.general_stats['headers']['featureCounts_percent'] = '<th class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="featureCounts: % Assigned reads">% Assigned</span></th>'
        config.general_stats['headers']['featureCounts'] = '<th class="chroma-col" data-chroma-scale="PuBu" data-chroma-min="0"><span data-toggle="tooltip" title="featureCounts: Assigned reads (millions)">M Assigned</span></th>'
        for sn, data in self.featurecounts_data.items():
            config.general_stats['rows'][sn]['featureCounts_percent'] = '<td class="text-right">{:.1f}%</td>'.format((data['Assigned']/data['Total'])*100)
            config.general_stats['rows'][sn]['featureCounts'] = '<td class="text-right">{:.1f}</td>'.format(data['Assigned']/1000000)


    def featureCounts_chart (self):
        """ Make the featureCounts assignment rates plot """
        
        # Specify the order of the different possible categories
        keys = ['Assigned', 'Unassigned_Ambiguity', 'Unassigned_MultiMapping', 'Unassigned_NoFeatures',
        'Unassigned_Unmapped', 'Unassigned_MappingQuality', 'Unassigned_FragementLength', 'Unassigned_Chimera',
        'Unassigned_Secondary', 'Unassigned_Nonjunction', 'Unassigned_Duplicate']
        
        # Config for the plot
        config = {
            'id': 'featureCounts_assignment_plot',
            'title': 'featureCounts Assignments',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.featurecounts_data, keys, config)
