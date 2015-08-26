#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
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
        self.cutadapt_data = dict()
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                if os.path.getsize(os.path.join(root,fn)) < 10000:
                    try:
                        with open (os.path.join(root,fn), "r") as f:
                            s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                            self.parse_cutadapt_logs(s)
                    except ValueError:
                        logging.debug("Couldn't read file when looking for cutadapt output: {}".format(fn))

        if len(self.cutadapt_data) == 0:
            logging.debug("Could not find any Cutadapt reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Cutadapt reports".format(len(self.cutadapt_data)))

        # Write parsed report data to a file
        # Only the summary stats - skip the length data (t_lengths)
        with open (os.path.join(self.output_dir, 'report_data', 'multiqc_cutadapt.txt'), "w") as f:
            print( self.dict_to_csv( { k: { j: x for j, x in v.iteritems() if j != 't_lengths'} for k, v in self.cutadapt_data.iteritems() } ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.cutadapt_general_stats_table(report)

        # Trimming Length Profiles
        # Only one section, so add to the intro
        self.intro += self.cutadapt_length_trimmed_plot()


    def parse_cutadapt_logs(self, s):
        i = s.find('Input filename', 0)
        while i >= 0:
            s.seek(i)

            fn_search = re.search(r"^Input filename:\s+(.+)$", s.readline())
            while fn_search is None:
                l = s.readline()
                if not l: return # end of file
                fn_search = re.search(r"^Input filename:\s+(.+)$", l)
            s_name = fn_search.group(1)
            s_name = s_name.split(".txt",1)[0]
            s_name = s_name.split("_trimming_report",1)[0]
            s_name = s_name.split(".gz",1)[0]
            s_name = s_name.split(".fastq",1)[0]
            s_name = s_name.split(".fq",1)[0]

            self.cutadapt_data[s_name] = dict()
            regexes = {
                'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
                'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'r_processed': "Total reads processed:\s*([\d,]+)",
                'r_with_adapters': "Reads with adapters:\s*([\d,]+)"
            }

            l = s.readline()
            while l:
                for k, r in regexes.iteritems():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(',', ''))

                if l.find('length') > -1 and l.find('max.err') > -1:
                    l = s.readline()
                    r_seqs = re.search(r"^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                    self.cutadapt_data[s_name]['t_lengths'] = {}
                    while r_seqs:
                        a_len = int(r_seqs.group(1))
                        self.cutadapt_data[s_name]['t_lengths'][a_len] = {
                            'count': int(r_seqs.group(2)),
                            'expect': float(r_seqs.group(3))
                        }
                        if float(r_seqs.group(3)) > 0:
                            self.cutadapt_data[s_name]['t_lengths'][a_len]['obs_exp'] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                        else:
                            # Cheating, I know. Infinity is difficult to plot.
                            self.cutadapt_data[s_name]['t_lengths'][a_len]['obs_exp'] = float(r_seqs.group(2))

                        l = s.readline() # Push on to next line for regex
                        r_seqs = re.search(r"^(\d+)\s+(\d+)\s+([\d\.]+)", l)

                if len(self.cutadapt_data[s_name]) == len(regexes) + 1: break # Found everything we need
                if l.find('Input filename') > -1: break # Didn't find everything, but come across another sample

                l = s.readline() # Next line for while loop

            if 'bp_processed' in self.cutadapt_data[s_name] and 'bp_written' in self.cutadapt_data[s_name]:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(self.cutadapt_data[s_name]['bp_processed'] - self.cutadapt_data[s_name]['bp_written']) / self.cutadapt_data[s_name]['bp_processed']) * 100

            # Look for the next cutadapt output in this file
            i = s.find('This is cutadapt', i + 1)


    def cutadapt_general_stats_table(self, report):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bp_trimmed'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="% Total Base Pairs trimmed by Cutadapt">Trimmed</span></th>'
        for samp, vals in self.cutadapt_data.iteritems():
            report['general_stats']['rows'][samp]['bp_trimmed'] = '<td class="text-right">{:.1f}%</td>'.format(vals['percent_trimmed'])

    def cutadapt_length_trimmed_plot (self):

        counts = list()
        obsexp = list()
        for s in sorted(self.cutadapt_data):
            counts_pairs = list()
            obsexp_pairs = list()
            for l, p in iter(sorted(self.cutadapt_data[s]['t_lengths'].iteritems())):
                counts_pairs.append([l, p['count']])
                obsexp_pairs.append([l, p['obs_exp']])
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
