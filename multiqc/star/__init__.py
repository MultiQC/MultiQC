#!/usr/bin/env python

""" MultiQC module to parse output from STAR """

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
        self.name = "STAR"
        self.anchor = "star"
        self.intro = '<p><a href="https://github.com/alexdobin/STAR" target="_blank">STAR</a> \
            is an ultrafast universal RNA-seq aligner.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any STAR reports
        self.star_data = defaultdict(lambda:dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith('Log.final.out'):
                    with open (os.path.join(root,fn), "r") as f:
                        parsed_data = self.parse_star_report(f.read())
                        if parsed_data is not None:
                            self.star_data[fn[:-13]] = parsed_data

        if len(self.star_data) == 0:
            logging.debug("Could not find any STAR reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} STAR reports".format(len(self.star_data)))

        # Write parsed report data to a file
        with open (os.path.join(self.output_dir, 'report_data', 'multiqc_star.txt'), "w") as f:
            print( self.dict_to_csv( self.star_data ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.star_stats_table(report)

        # Alignment bar plot
        # Only one section, so add to the intro
        self.intro += self.star_alignment_chart()


    def parse_star_report (self, raw_data):
        """ Parse the final STAR log file. """

        regexes = {
            'total_reads':                  r"Number of input reads \|\s+(\d+)",
            'avg_input_read_length':        r"Average input read length \|\s+([\d\.]+)",
            'uniquely_mapped':              r"Uniquely mapped reads number \|\s+(\d+)",
            'uniquely_mapped_percent':      r"Uniquely mapped reads % \|\s+([\d\.]+)",
            'avg_mapped_read_length':       r"Average mapped length \|\s+([\d\.]+)",
            'num_splices':                  r"Number of splices: Total \|\s+(\d+)",
            'num_annotated_splices':        r"Number of splices: Annotated \(sjdb\) \|\s+(\d+)",
            'num_GTAG_splices':             r"Number of splices: GT/AG \|\s+(\d+)",
            'num_GCAG_splices':             r"Number of splices: GC/AG \|\s+(\d+)",
            'num_ATAC_splices':             r"Number of splices: AT/AC \|\s+(\d+)",
            'num_noncanonical_splices':     r"Number of splices: Non-canonical \|\s+(\d+)",
            'mismatch_rate':                r"Mismatch rate per base, % \|\s+([\d\.]+)",
            'deletion_rate':                r"Deletion rate per base \|\s+([\d\.]+)",
            'deletion_length':              r"Deletion average length \|\s+([\d\.]+)",
            'insertion_rate':               r"Insertion rate per base \|\s+([\d\.]+)",
            'insertion_length':             r"Insertion average length \|\s+([\d\.]+)",
            'multimapped':                  r"Number of reads mapped to multiple loci \|\s+(\d+)",
            'multimapped_percent':          r"% of reads mapped to multiple loci \|\s+([\d\.]+)",
            'multimapped_toomany':          r"Number of reads mapped to too many loci \|\s+(\d+)",
            'multimapped_toomany_percent':  r"% of reads mapped to too many loci \|\s+([\d\.]+)",
            'unmapped_mismatches_percent':  r"% of reads unmapped: too many mismatches \|\s+([\d\.]+)",
            'unmapped_tooshort_percent':    r"% of reads unmapped: too short \|\s+([\d\.]+)",
            'unmapped_other_percent':       r"% of reads unmapped: other \|\s+([\d\.]+)",
        }
        parsed_data = {}
        for k, r in regexes.iteritems():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
        # Figure out the numbers for unmapped as for some reason only the percentages are given
        try:
            total_mapped = parsed_data['uniquely_mapped'] + parsed_data['multimapped'] + parsed_data['multimapped_toomany']
            unmapped_count = parsed_data['total_reads'] - total_mapped
            total_unmapped_percent = parsed_data['unmapped_mismatches_percent'] + parsed_data['unmapped_tooshort_percent'] + parsed_data['unmapped_other_percent']
            parsed_data['unmapped_mismatches'] = int(round(unmapped_count * (parsed_data['unmapped_mismatches_percent'] / total_unmapped_percent), 0))
            parsed_data['unmapped_tooshort'] = int(round(unmapped_count * (parsed_data['unmapped_tooshort_percent'] / total_unmapped_percent), 0))
            parsed_data['unmapped_other'] = int(round(unmapped_count * (parsed_data['unmapped_other_percent'] / total_unmapped_percent), 0))
        except KeyError:
            pass

        if len(parsed_data) == 0: return None
        return parsed_data


    def star_stats_table(self, report):
        """ Take the parsed stats from the STAR report and add them to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['uniquely_mapped_percent'] = '<th class="chroma-col" data-chroma-scale="YlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="STAR: % Uniquely mapped reads">%&nbsp;Mapped</span></th>'
        report['general_stats']['headers']['uniquely_mapped'] = '<th class="chroma-col" data-chroma-scale="PuRd" data-chroma-min="0"><span data-toggle="tooltip" title="STAR: Uniquely mapped reads (millions)">M&nbsp;Mapped</span></th>'
        for sn, data in self.star_data.iteritems():
            report['general_stats']['rows'][sn]['uniquely_mapped_percent'] = '<td class="text-right">{:.1f}%</td>'.format(data['uniquely_mapped_percent'])
            report['general_stats']['rows'][sn]['uniquely_mapped'] = '<td class="text-right">{:.1f}</td>'.format(data['uniquely_mapped']/1000000)

    def star_alignment_chart (self):
        """ Make the HighCharts HTML to plot the alignment rates """

        cats = sorted(self.star_data.keys())
        data = list()
        keys = OrderedDict()
        keys['unmapped_other'] =      'Unmapped: other'
        keys['unmapped_tooshort'] =   'Unmapped: too short'
        keys['unmapped_mismatches'] = 'Unmapped: too many mismatchess'
        keys['multimapped_toomany'] = 'Mapped to too many loci'
        keys['multimapped'] =         'Mapped to multiple loci'
        keys['uniquely_mapped'] =     'Uniquely mapped'
        colours = {
            'unmapped_other':      '#7f0000',
            'unmapped_tooshort':   '#b1084c',
            'unmapped_mismatches': '#e63491',
            'multimapped_toomany': '#f7a35c',
            'multimapped':         '#7cb5ec',
            'uniquely_mapped':     '#437bb1',
        }

        for k, name in keys.iteritems():
            thisdata = list()
            for sn in cats:
                thisdata.append(self.star_data[sn][k])
            if max(thisdata) > 0:
                data.append({
                    'name': name,
                    'color': colours[k],
                    'data': thisdata
                })

        return '<div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#star_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_percent" data-target="#star_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="star_alignment_plot" class="fastqc-overlay-plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            star_alignment_cats = {};\n\
            star_alignment_data = {};\n\
            var star_alignment_pconfig = {{ \n\
                "title": "STAR Alignment Scores",\n\
                "ylab": "# Reads",\n\
                "ymin": 0,\n\
                "stacking": "normal" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#star_alignment_plot", star_alignment_cats, star_alignment_data, star_alignment_pconfig); \
            }}); \
        </script>'.format(json.dumps(cats), json.dumps(data));
