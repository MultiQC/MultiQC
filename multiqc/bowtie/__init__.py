#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie """

from __future__ import print_function
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
        self.anchor = "bowtie"
        self.intro = '<p><a href="http://bowtie-bio.sourceforge.net/" target="_blank">Bowtie</a> \
            is an ultrafast, memory-efficient short read aligner.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bowtie reports
        self.bowtie_data = dict()
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                if os.path.getsize(os.path.join(root,fn)) < 10000:
                    try:
                        with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                            s = f.read()
                            parsed_data = self.parse_bowtie_logs(s)
                            if parsed_data is not None:
                                self.bowtie_data[fn] = parsed_data
                    except ValueError:
                        logging.debug("Couldn't read file when looking for bowtie output: {}".format(fn))

        if len(self.bowtie_data) == 0:
            logging.debug("Could not find any Bowtie reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Bowtie reports".format(len(self.bowtie_data)))

        # Write parsed report data to a file
        # Only the summary stats - skip the length data (t_lengths)
        with open (os.path.join(self.output_dir, 'report_data', 'multiqc_bowtie.txt'), "w") as f:
            print( self.dict_to_csv( { k: { j: x for j, x in v.items() if j != 't_lengths'} for k, v in self.bowtie_data.items() } ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie_general_stats_table(report)

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.bowtie_alignment_plot()


    def parse_bowtie_logs(self, s):
        i = s.find('# reads processed:', 0)
        parsed_data = {}
        if i >= 0:
            regexes = {
                'reads_processed': r"# reads processed:\s+(\d+)",
                'reads_aligned': r"# reads with at least one reported alignment:\s+(\d+)",
                'reads_aligned_percentage': r"# reads with at least one reported alignment:\s+\d+\s+\(([\d\.]+)%\)",
                'not_aligned': r"# reads that failed to align:\s+(\d+)",
                'not_aligned_percentage': r"# reads that failed to align:\s+\d+\s+\(([\d\.]+)%\)",
                'mulimapped': r"# reads with alignments suppressed due to -m:\s+(\d+)",
                'mulimapped_percentage': r"# reads with alignments suppressed due to -m:\s+\d+\s+\(([\d\.]+)%\)"
            }

            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = int(match.group(1).replace(',', ''))
        if len(parsed_data) == 0: return None
        return parsed_data


    def bowtie_general_stats_table(self, report):
        """ Take the parsed stats from the Bowtie report and add it to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bowtie_aligned'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bowtie: % reads with at least one reported alignment">%&nbsp;Aligned</span></th>'
        for samp, vals in self.bowtie_data.items():
            report['general_stats']['rows'][samp]['bowtie_aligned'] = '<td class="text-right">{:.1f}%</td>'.format(vals['reads_aligned_percentage'])

    def bowtie_alignment_plot (self):
        return ''
        # counts = list()
        # obsexp = list()
        # for s in sorted(self.bowtie_data):
        #     counts_pairs = list()
        #     obsexp_pairs = list()
        #     for l, p in iter(sorted(self.bowtie_data[s]['t_lengths'].items())):
        #         counts_pairs.append([l, p['count']])
        #         obsexp_pairs.append([l, p['obs_exp']])
        #     counts.append({
        #         'name': s,
        #         'data': counts_pairs
        #     })
        #     obsexp.append({
        #         'name': s,
        #         'data': obsexp_pairs
        #     })
        #
        # html = '<p>This plot shows the number of reads with certain lengths of adapter trimmed. \n\
        # Obs/Exp shows the raw counts divided by the number expected due to sequencing errors. A defined peak \n\
        # may be related to adapter length. See the \n\
        # <a href="http://bowtie.readthedocs.org/en/latest/guide.html#how-to-read-the-report" target="_blank">bowtie documentation</a> \n\
        # for more information on how these numbers are generated.</p> \n\
        # <div class="btn-group switch_group"> \n\
		# 	<button class="btn btn-default btn-sm active" data-action="set_data" data-ylab="Obs / Expected" data-newdata="bowtie_length_obsexp" data-target="#bowtie_length_plot">Obs/Exp</button> \n\
		# 	<button class="btn btn-default btn-sm" data-action="set_data" data-ylab="Count" data-newdata="bowtie_length_counts" data-target="#bowtie_length_plot">Counts</button> \n\
		# </div> \n\
        # <div id="bowtie_length_plot" class="hc-plot"></div> \n\
        # <script type="text/javascript"> \n\
        #     bowtie_length_counts = {};\n\
        #     bowtie_length_obsexp = {};\n\
        #     var bowtie_length_pconfig = {{ \n\
        #         "title": "Lengths Trimmed",\n\
        #         "ylab": "Obs / Expected",\n\
        #         "xlab": "Length Trimmed (bp)",\n\
        #         "ymin": 0,\n\
        #         "tt_label": "<b>{{point.x}} bp trimmed</b>",\n\
        #         "use_legend": false,\n\
        #     }}; \n\
        #     $(function () {{ \
        #         plot_xy_line_graph("#bowtie_length_plot", bowtie_length_obsexp, bowtie_length_pconfig); \n\
        #     }}); \
        # </script>'.format(json.dumps(counts), json.dumps(obsexp));
        #
        # return html
