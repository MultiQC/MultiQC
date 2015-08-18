#!/usr/bin/env python

""" MultiQC module to parse output from Bismark """

from collections import defaultdict
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
        self.name = "Bismark"
        self.anchor = "bismark"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/bismark/" target="_blank">Bismark</a> \
            is a tool to map bisulfite converted sequence reads and determine cytosine methylation states.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bismark reports
        self.bismark_raw_data = defaultdict(lambda:dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith('_PE_report.txt') or fn.endswith('_SE_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        fn_search = re.search("Bismark report for: (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            self.bismark_raw_data[s_name]['alignment'] = r_data
                        else:
                            logging.warn("Found bismark alignment report, but couldn't recognise contents: {}".format(fn))

                if fn.endswith('.deduplication_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        fn_search = re.search("Total number of alignments analysed in (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            self.bismark_raw_data[s_name]['dedup'] = r_data
                        else:
                            logging.warn("Found bismark deduplication report, but couldn't recognise contents: {}".format(fn))

                if fn.endswith('_splitting_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        s_name = r_data.splitlines()[0]
                        s_name = s_name.split(".gz",1)[0]
                        s_name = s_name.split(".fastq",1)[0]
                        s_name = s_name.split(".fq",1)[0]
                        s_name = s_name.split("_val_1",1)[0]
                        self.bismark_raw_data[s_name]['methextract'] = r_data

        if len(self.bismark_raw_data) == 0:
            logging.debug("Could not find any Bismark reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Bismark reports".format(len(self.bismark_raw_data)))

        # Parse the raw reports
        self.parse_bismark_reports()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bismark_stats_table(report)

        self.sections = list()

    def parse_bismark_reports (self):
        """ Search the three types of Bismark report files for
        numbers needed later in the module. """

        regexes = {
            'alignment': {
                'total_reads': r"^Sequence(?:s| pairs) analysed in total:\s+(\d+)$",
                'no_alignments': r"^Sequence(?:s| pairs) with no alignments under any condition:\s+(\d+)$",
                'ambig_reads': r"^Sequence(?:s| pairs) did not map uniquely:\s+(\d+)$",
                'discarded_reads': r"^Sequence(?:s| pairs) which were discarded because genomic sequence could not be extracted:\s+(\d+)$"
            },
            'dedup': {
                'aligned_reads': r"^Total number of alignments analysed in .+:\s+(\d+)$",
                'dup_reads': r"^Total number duplicated alignments removed:\s+(\d+)",
                'dup_reads_percent': r"^Total number duplicated alignments removed:\s+\d+\s+\(([\d\.]+)%\)",
                'dedup_reads': r"^Total count of deduplicated leftover sequences:\s+(\d+)",
                'dedup_reads_percent': r"^Total count of deduplicated leftover sequences:\s+\d+\s+\(([\d\.]+)% of total\)"
            },
            'methextract': {
                'total_c': r"^Total number of C's analysed:\s+(\d+)$",
                'meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
                'meth_cph': r"^Total methylated C's in CHG context:\s+(\d+)",
                'meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
                'unmeth_cpg': r"^Total C to T conversions in CpG context:\s+(\d+)",
                'unmeth_cph': r"^Total C to T conversions in CHG context:\s+(\d+)",
                'unmeth_chh': r"^Total C to T conversions in CHH context:\s+(\d+)"
            }
        }
        for sn, data in self.bismark_raw_data.iteritems():
            for report_type in regexes.keys():
                for k, r in regexes[report_type].iteritems():
                    r_search = re.search(r, data[report_type], re.MULTILINE)
                    if r_search:
                        self.bismark_raw_data[sn][k] = float(r_search.group(1))

    def bismark_stats_table(self, report):
        """ Take the parsed stats from the Bismark reports and add them to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bismark_dedup_reads'] = '<th class="chroma-col" data-chroma-scale="Greens" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Deduplicated Alignments (millions)">M Unique</span></th>'
        report['general_stats']['headers']['bismark_dedup_reads_percent'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Percent Duplicated Alignments">% Dups</span></th>'
        report['general_stats']['headers']['bismark_aligned'] = '<th class="chroma-col" data-chroma-scale="PuRd" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total Aligned Sequences (millions)">M Aligned</span></th>'

        for sn, data in self.bismark_raw_data.iteritems():
            report['general_stats']['rows'][sn]['bismark_dedup_reads'] = '<td class="text-right">{:.1f}</td>'.format(data['dedup_reads']/1000000)
            report['general_stats']['rows'][sn]['bismark_dedup_reads_percent'] = '<td class="text-right">{}%</td>'.format(data['dup_reads_percent'])
            report['general_stats']['rows'][sn]['bismark_aligned'] = '<td class="text-right">{:.1f}</td>'.format(data['aligned_reads']/1000000)
