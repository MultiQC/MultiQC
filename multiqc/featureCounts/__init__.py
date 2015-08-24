#!/usr/bin/env python

""" MultiQC module to parse output from featureCounts """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import json
import logging
import mmap
import os
import re

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "featureCounts"
        self.anchor = "featurecounts"
        self.intro = '<p><a href="http://bioinf.wehi.edu.au/featureCounts/" target="_blank">Subread featureCounts</a> \
             is a highly efficient general-purpose read summarization program that counts mapped reads for genomic \
             features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any featureCounts reports
        self.featurecounts_data = defaultdict(lambda:dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith('_counts.txt.summary'):
                    with open (os.path.join(root,fn), "r") as f:
                        s_name = fn[:-19]
                        s_name = s_name.split(".bam",1)[0]
                        s_name = s_name.split(".sam",1)[0]
                        s_name = s_name.split("_star_aligned",1)[0]
                        s_name = s_name.split(".gz",1)[0]
                        s_name = s_name.split(".fastq",1)[0]
                        s_name = s_name.split(".fq",1)[0]
                        parsed_data = self.parse_featurecounts_report(f.read())
                        if parsed_data is not None:
                            self.featurecounts_data[s_name] = parsed_data

        if len(self.featurecounts_data) == 0:
            logging.debug("Could not find any featureCounts reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} featureCounts reports".format(len(self.featurecounts_data)))

        # Write parsed report data to a file
        with open (os.path.join(self.output_dir, 'report_data', 'multiqc_featureCounts.txt'), "w") as f:
            print( self.dict_to_csv( self.featurecounts_data ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.featurecounts_stats_table(report)

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
        for k, r in regexes.iteritems():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
                total_count += float(r_search.group(1))
        parsed_data['Total'] = total_count
        if len(parsed_data) == 0: return None
        return parsed_data


    def featurecounts_stats_table(self, report):
        """ Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['featureCounts_percent'] = '<th class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="featureCounts: % Assigned reads">%&nbsp;Assigned</span></th>'
        report['general_stats']['headers']['featureCounts'] = '<th class="chroma-col" data-chroma-scale="PuBu" data-chroma-min="0"><span data-toggle="tooltip" title="featureCounts: Assigned reads (millions)">M&nbsp;Assigned</span></th>'
        for sn, data in self.featurecounts_data.iteritems():
            report['general_stats']['rows'][sn]['featureCounts_percent'] = '<td class="text-right">{:.1f}%</td>'.format((data['Assigned']/data['Total'])*100)
            report['general_stats']['rows'][sn]['featureCounts'] = '<td class="text-right">{:.1f}</td>'.format(data['Assigned']/1000000)


    def featureCounts_chart (self):
        """ Make the HighCharts HTML to plot the assignment rates """

        cats = sorted(self.featurecounts_data.keys())
        data = list()
        keys = ['Assigned', 'Unassigned_Ambiguity', 'Unassigned_MultiMapping', 'Unassigned_NoFeatures',
        'Unassigned_Unmapped', 'Unassigned_MappingQuality', 'Unassigned_FragementLength', 'Unassigned_Chimera',
        'Unassigned_Secondary', 'Unassigned_Nonjunction', 'Unassigned_Duplicate']

        for k in keys:
            thisdata = list()
            name = k.replace('_', ' ')
            for sn in cats:
                thisdata.append(self.featurecounts_data[sn][k])
            if max(thisdata) > 0:
                data.append({
                    'name': name,
                    'data': thisdata
                })

        return '<div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#featurecounts_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_percent" data-target="#featurecounts_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="featurecounts_alignment_plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            featurecounts_alignment_cats = {};\n\
            featurecounts_alignment_data = {};\n\
            var featurecounts_alignment_pconfig = {{ \n\
                "title": "featurecounts Alignment Scores",\n\
                "ylab": "# Reads",\n\
                "ymin": 0,\n\
                "stacking": "normal" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#featurecounts_alignment_plot", featurecounts_alignment_cats, featurecounts_alignment_data, featurecounts_alignment_pconfig); \
            }}); \
        </script>'.format(json.dumps(cats), json.dumps(data));
