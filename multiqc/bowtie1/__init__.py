#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 1 """

from __future__ import print_function
from collections import OrderedDict
import io
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
        self.name = "Bowtie"
        self.anchor = "bowtie1"
        self.intro = '<p><a href="http://bowtie-bio.sourceforge.net/" target="_blank">Bowtie 1</a> \
            is an ultrafast, memory-efficient short read aligner.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bowtie reports
        self.bowtie_data = dict()
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                if os.path.getsize(os.path.join(root,fn)) < 200000:
                    try:
                        with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                            s = f.read()
                            parsed_data = self.parse_bowtie_logs(s)
                            if parsed_data is not None:
                                s_name = self.clean_s_name(fn, root, prepend_dirs=report['prepend_dirs'])
                                if s_name in self.bowtie_data:
                                    logging.warn("Bowtie 1: Duplicate sample name found! Overwriting: {}".format(s_name))
                                self.bowtie_data[s_name] = parsed_data
                    except ValueError:
                        logging.debug("Couldn't read file when looking for Bowtie 1 output: {}".format(fn))

        if len(self.bowtie_data) == 0:
            logging.debug("Could not find any Bowtie 1 reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Bowtie 1 reports".format(len(self.bowtie_data)))

        # Write parsed report data to a file
        with io.open (os.path.join(self.output_dir, 'report_data', 'multiqc_bowtie1.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( { k: { j: x for j, x in v.items() if j != 't_lengths'} for k, v in self.bowtie_data.items() } ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie_general_stats_table(report)

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.bowtie_alignment_plot()


    def parse_bowtie_logs(self, s):
        # Check that this isn't actually Bismark using bowtie
        if s.find('Using bowtie 1 for aligning with bismark.', 0) >= 0: return None
        i = s.find('# reads processed:', 0)
        parsed_data = {}
        if i >= 0:
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
        if len(parsed_data) == 0: return None
        return parsed_data


    def bowtie_general_stats_table(self, report):
        """ Take the parsed stats from the Bowtie report and add it to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bowtie_aligned'] = '<th class="chroma-col" data-chroma-scale="OrRd-rev" data-chroma-max="100" data-chroma-min="20"><span data-toggle="tooltip" title="Bowtie 1: % reads with at least one reported alignment">%&nbsp;Aligned</span></th>'
        for samp, vals in self.bowtie_data.items():
            report['general_stats']['rows'][samp]['bowtie_aligned'] = '<td class="text-right">{:.1f}%</td>'.format(vals['reads_aligned_percentage'])

    def bowtie_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """

        cats = sorted(self.bowtie_data.keys())
        data = list()
        keys = OrderedDict()
        keys['reads_aligned'] = 'Aligned'
        keys['multimapped'] =    'Multimapped'
        keys['not_aligned'] =   'Not aligned'

        colours = {
            'reads_aligned':  '#8bbc21',
            'multimapped':    '#2f7ed8',
            'not_aligned':    '#0d233a',
        }

        for k, name in keys.items():
            thisdata = list()
            for sn in cats:
                thisdata.append(self.bowtie_data[sn].get(k, 0))
            if max(thisdata) > 0:
                data.append({
                    'name': name,
                    'color': colours[k],
                    'data': thisdata
                })

        return '<div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#bowtie1_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_percent" data-target="#bowtie1_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="bowtie1_alignment_plot" class="hc-plot"></div> \n\
        <script type="text/javascript"> \n\
            bowtie1_alignment_cats = {};\n\
            bowtie1_alignment_data = {};\n\
            var bowtie1_alignment_pconfig = {{ \n\
                "title": "Bowtie 1 Alignment Scores",\n\
                "ylab": "# Reads",\n\
                "ymin": 0,\n\
                "stacking": "normal" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bowtie1_alignment_plot", bowtie1_alignment_cats, bowtie1_alignment_data, bowtie1_alignment_pconfig); \
            }}); \
        </script>'.format(json.dumps(cats), json.dumps(data));
