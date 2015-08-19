#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

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
        self.name = "Cutadapt"
        self.anchor = "cutadapt"
        self.intro = '<p><a href="https://code.google.com/p/cutadapt/" target="_blank">Cutadapt</a> \
            is a tool to find and remove adapter sequences, primers, poly-A tails and other types \
            of unwanted sequence from your high-throughput sequencing reads.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Cutadapt reports
        cutadapt_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith('.txt') and os.path.getsize(os.path.join(root,fn)) < 50000:
                    with open (os.path.join(root,fn), "r") as f:
                        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                        if s.find('This is cutadapt') != -1:
                            fn_search = re.search(r"^Input filename:\s+(.+)$", s, re.MULTILINE)
                            if fn_search:
                                s_name = fn_search.group(1)
                            else:
                                s_name = fn
                            s_name = s_name.split(".txt",1)[0]
                            s_name = s_name.split("_trimming_report",1)[0]
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]

                            cutadapt_raw_data[s_name] = f.read()

        if len(cutadapt_raw_data) == 0:
            logging.debug("Could not find any Cutadapt reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Cutadapt reports".format(len(cutadapt_raw_data)))

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        parsed_stats = self.cutadapt_general_stats(cutadapt_raw_data)
        self.cutadapt_general_stats_table(parsed_stats, report)

        # Trimming Length Profiles
        # Only one section, so add to the intro
        length_trimmed = self.cutadapt_length_trimmed(cutadapt_raw_data)
        self.intro += self.cutadapt_length_trimmed_plot(length_trimmed)


    def cutadapt_general_stats(self, cutadapt_raw_data):
        """ Parse the single-digit stats for each sample from the Cutadapt report. """
        parsed_stats = {}
        for s, data in cutadapt_raw_data.iteritems():
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

    def cutadapt_general_stats_table(self, parsed_stats, report):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bp_trimmed'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="% Total Base Pairs trimmed by Cutadapt">Trimmed</span></th>'
        for samp, vals in parsed_stats.iteritems():
            report['general_stats']['rows'][samp]['bp_trimmed'] = '<td class="text-right">{:.1f}%</td>'.format(vals['percent_trimmed'])

    def cutadapt_length_trimmed(self, cutadapt_raw_data):
        """ Parse the counts of adapter lengths that have been trimmed """

        parsed_data = {}
        for s, data in cutadapt_raw_data.iteritems():
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
                        if float(r_seqs.group(3)) > 0:
                            parsed_data[s][a_len]['obs_exp'] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                        else:
                            parsed_data[s][a_len]['obs_exp'] = float(r_seqs.group(2)) # Cheating, I know. Infinity is difficult to plot.
                    else:
                        break
        return parsed_data

    def cutadapt_length_trimmed_plot (self, parsed_data):

        counts = list()
        obsexp = list()
        for s in sorted(parsed_data):
            counts_pairs = list()
            obsexp_pairs = list()
            for l, p in iter(sorted(parsed_data[s].iteritems())):
                counts_pairs.append([l, parsed_data[s][l]['count']])
                obsexp_pairs.append([l, parsed_data[s][l]['obs_exp']])
            counts.append({
                'name': s,
                'data': counts_pairs
            })
            obsexp.append({
                'name': s,
                'data': obsexp_pairs
            })

        html = '<p>This plot shows the number of reads with certain lengths of adapter trimmed. \n\
        Obs/Exp shows the raw counts divided by the number expected due to sequencing errors. A defined peak \n\
        may be related to adapter length. See the \n\
        <a href="http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report" target="_blank">cutadapt documentation</a> \n\
        for more information on how these numbers are generated.</p> \n\
        <div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_data" data-ylab="Obs / Expected" data-newdata="cutadapt_length_obsexp" data-target="#cutadapt_length_plot">Obs/Exp</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_data" data-ylab="Count" data-newdata="cutadapt_length_counts" data-target="#cutadapt_length_plot">Counts</button> \n\
		</div> \n\
        <div id="cutadapt_length_plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            cutadapt_length_counts = {};\n\
            cutadapt_length_obsexp = {};\n\
            var cutadapt_length_pconfig = {{ \n\
                "title": "Lengths Trimmed",\n\
                "ylab": "Obs / Expected",\n\
                "xlab": "Length Trimmed (bp)",\n\
                "ymin": 0,\n\
                "tt_label": "<b>{{point.x}} bp trimmed</b>",\n\
                "use_legend": false,\n\
            }}; \n\
            $(function () {{ \
                plot_xy_line_graph("#cutadapt_length_plot", cutadapt_length_obsexp, cutadapt_length_pconfig); \n\
            }}); \
        </script>'.format(json.dumps(counts), json.dumps(obsexp));

        return html
