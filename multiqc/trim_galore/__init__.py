#!/usr/bin/env python

""" MultiQC module to parse output from Trim Galore!
NB: This is actually cutadapt output that we're parsing, but there's
no fixed filename that we can search for with cutadapt.
"""

import json
import logging
import os
import re

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "Trim Galore!"
        self.anchor = "trim_galore"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">Trim Galore!</a> \
            is a wrapper tool around <a href="https://code.google.com/p/cutadapt/" target="_blank">cutadapt</a> \
            to consistently apply quality and adapter trimming to FastQ files.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Trim Galore! reports
        tgalore_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith("_trimming_report.txt"):
                    s_name = fn[:-20]
                    if s_name[-3:] == ".gz":
                        s_name = s_name[:-3]
                    if s_name[-6:] == ".fastq":
                        s_name = s_name[:-6]
                    if s_name[-3:] == ".fq":
                        s_name = s_name[:-3]
                    with open (os.path.join(root,fn), "r") as f:
                        tgalore_raw_data[s_name] = f.read()

        if len(tgalore_raw_data) == 0:
            logging.debug("Could not find any Trim Galore! reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Trim Galore! reports".format(len(tgalore_raw_data)))

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        parsed_stats = self.tgalore_basic_stats(tgalore_raw_data)
        self.tgalore_basic_stats_table(parsed_stats, report)

        # Section 1 - Trimming Length Profiles
        length_trimmed = self.trimgalore_length_trimmed(tgalore_raw_data)
        self.sections.append({
            'name': 'Trimming Length Profiles',
            'anchor': 'tgalore-lengths',
            'content': self.trimgalore_length_trimmed_plot(length_trimmed)
        })


    def tgalore_basic_stats(self, tgalore_raw_data):
        """ Parse the single-digit stats for each sample from the Trim Galore!
        report. """
        parsed_stats = {}
        for s, data in tgalore_raw_data.iteritems():
            parsed_stats[s] = {}

            bp_processed = re.search("Total basepairs processed:\s*([\d,]+) bp", data)
            if bp_processed:
                parsed_stats[s]['bp_processed'] = int(bp_processed.group(1).replace(',', ''))

            bp_written = re.search("Total written \(filtered\):\s*([\d,]+) bp", data)
            if bp_written:
                parsed_stats[s]['bp_written'] = int(bp_written.group(1).replace(',', ''))

            if bp_processed and bp_written:
                parsed_stats[s]['percent_trimmed'] = (float(parsed_stats[s]['bp_processed'] - parsed_stats[s]['bp_written']) / parsed_stats[s]['bp_processed']) * 100

            quality_trimmed = re.search("Quality-trimmed:\s*([\d,]+) bp", data)
            if quality_trimmed:
                parsed_stats[s]['quality_trimmed'] = int(quality_trimmed.group(1).replace(',', ''))

            r_processed = re.search("Total reads processed:\s*([\d,]+)", data)
            if r_processed:
                parsed_stats[s]['r_processed'] = int(r_processed.group(1).replace(',', ''))

            r_with_adapters = re.search("Reads with adapters:\s*([\d,]+)", data)
            if r_with_adapters:
                parsed_stats[s]['r_with_adapters'] = int(r_with_adapters.group(1).replace(',', ''))

        return parsed_stats

    def tgalore_basic_stats_table(self, parsed_stats, report):
        """ Take the parsed stats from the Trim Galore! report and add it to the
        basic stats table at the top of the report """

        report['basic_stats']['headers']['bp_trimmed'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="% Total Base Pairs trimmed by Trim Galore!">Trimmed</span></th>'
        for samp, vals in parsed_stats.iteritems():
            report['basic_stats']['rows'][samp]['bp_trimmed'] = '<td class="text-right">{:.1f}%</td>'.format(vals['percent_trimmed'])

    def trimgalore_length_trimmed(self, tgalore_raw_data):
        """ Parse the 'Per base sequence quality' data from fastqc_data.txt
        Returns a 2D dict, sample names as first keys, then a dict of lists
        containing base, mean, median, lower_quart, upper_quart, 10_percentile
        and 90_percentile. """

        parsed_data = {}
        for s, data in tgalore_raw_data.iteritems():
            parsed_data[s] = {}
            in_section = False
            for l in data.splitlines():
                if l == "length	count	expect	max.err	error counts":
                    in_section = True
                elif in_section:
                    r_seqs = re.search(r"^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                    if r_seqs:
                        a_len = int(r_seqs.group(1))
                        parsed_data[s][a_len] = {}
                        parsed_data[s][a_len]['count'] = int(r_seqs.group(2))
                        parsed_data[s][a_len]['expect'] = float(r_seqs.group(3))
                        parsed_data[s][a_len]['obs_exp'] = int(r_seqs.group(2)) / float(r_seqs.group(3))
                    else:
                        break
        return parsed_data

    def trimgalore_length_trimmed_plot (self, parsed_data):

        data = list()
        for s in sorted(parsed_data):
            pairs = list()
            for l, p in iter(sorted(parsed_data[s].iteritems())):
                pairs.append([l, parsed_data[s][l]['obs_exp']])
            data.append({
                'name': s,
                'data': pairs
            })

        html = '<p>This plot shows the number of reads with certain lengths of adapter trimmed. \n\
        These counts are divided by the number expected due to sequencing errors. See the \n\
        <a href="http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report" target="_blank">cutadapt documentation</a> \n\
        for more information on how these numbers are generated.</p> \n\
        <div id="tgalore_length_trimmed" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            tgalore_length_trimmed_data = {};\n\
            var tgalore_l_pconfig = {{ \n\
                "title": "Lengths Trimmed",\n\
                "ylab": "Obs / Expected",\n\
                "xlab": "Length Trimmed (bp)",\n\
                "ymin": 0,\n\
                "tt_label": "<b>{{point.x}}bp trimmed</b>",\n\
                "use_legend": false,\n\
            }}; \n\
            $(function () {{ \
                plot_xy_line_graph("#tgalore_length_trimmed", tgalore_length_trimmed_data, tgalore_l_pconfig); \
            }}); \
        </script>'.format(json.dumps(data));

        return html
