#!/usr/bin/env python

""" MultiQC module to parse output from Bismark """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
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
        self.name = "Bismark"
        self.anchor = "bismark"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/bismark/" target="_blank">Bismark</a> \
            is a tool to map bisulfite converted sequence reads and determine cytosine methylation states.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bismark reports
        self.bismark_raw_data = defaultdict(lambda:dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                if fn.endswith('_PE_report.txt') or fn.endswith('_SE_report.txt'):
                    with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                        r_data = f.read()
                        fn_search = re.search("Bismark report for: (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            if report['prepend_dirs']:
                                s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
                            self.bismark_raw_data[s_name]['alignment'] = r_data
                        else:
                            logging.warn("Didn't recognise bismark alignment report contents: {}".format(fn))

                if fn.endswith('.deduplication_report.txt'):
                    with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                        r_data = f.read()
                        fn_search = re.search("Total number of alignments analysed in (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            if report['prepend_dirs']:
                                s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
                            self.bismark_raw_data[s_name]['dedup'] = r_data
                        else:
                            logging.warn("Didn't recognise bismark deduplication report contents: {}".format(fn))

                if fn.endswith('_splitting_report.txt'):
                    with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                        r_data = f.read()
                        s_name = r_data.splitlines()[0]
                        s_name = s_name.split(".gz",1)[0]
                        s_name = s_name.split(".fastq",1)[0]
                        s_name = s_name.split(".fq",1)[0]
                        s_name = s_name.split("_val_1",1)[0]
                        if report['prepend_dirs']:
                            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
                        self.bismark_raw_data[s_name]['methextract'] = r_data

        if len(self.bismark_raw_data) == 0:
            logging.debug("Could not find any Bismark reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Bismark reports".format(len(self.bismark_raw_data)))

        # Parse the raw reports
        self.bismark_data = defaultdict(lambda:dict())
        self.parse_bismark_reports()

        # Write parsed report data to a file
        with io.open (os.path.join(self.output_dir, 'report_data', 'multiqc_bismark.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( self.bismark_data ), file=f)

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bismark_stats_table(report)

        self.sections = list()

        # Section 1 - Column chart of alignment stats
        self.parse_alignment_chart_data()
        self.sections.append({
            'name': 'Alignment Rates',
            'anchor': 'bismark-alignment',
            'content': self.bismark_alignment_chart()
        })

        # Section 2 - Methylation percentages
        self.parse_methylation_chart_data()
        self.sections.append({
            'name': 'Cytosine Methylation',
            'anchor': 'bismark-methylation',
            'content': self.bismark_methlyation_chart()
        })

        # Section 3 - Strand Alignments
        self.parse_strand_chart_data()
        self.sections.append({
            'name': 'Strand Alignment',
            'anchor': 'bismark-strands',
            'content': self.bismark_strand_chart()
        })

    def parse_bismark_reports (self):
        """ Search the three types of Bismark report files for
        numbers needed later in the module. """

        regexes = {
            'alignment': {
                'total_reads': r"^Sequence(?:s| pairs) analysed in total:\s+(\d+)$",
                'aligned_reads': r"^Number of(?: paired-end)? alignments with a unique best hit(?: from the different alignments)?:\s+(\d+)$",
                'no_alignments': r"^Sequence(?:s| pairs) with no alignments under any condition:\s+(\d+)$",
                'ambig_reads': r"^Sequence(?:s| pairs) did not map uniquely:\s+(\d+)$",
                'discarded_reads': r"^Sequence(?:s| pairs) which were discarded because genomic sequence could not be extracted:\s+(\d+)$",
                'aln_total_c': r"^Total number of C's analysed:\s+(\d+)$",
                'aln_meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
                'aln_meth_chg': r"^Total methylated C's in CHG context:\s+(\d+)",
                'aln_meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
                'aln_unmeth_cpg': r"^Total unmethylated C's in CpG context:\s+(\d+)",
                'aln_unmeth_chg': r"^Total unmethylated C's in CHG context:\s+(\d+)",
                'aln_unmeth_chh': r"^Total unmethylated C's in CHH context:\s+(\d+)",
                'aln_percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
                'aln_percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
                'aln_percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%",
                'aln_strand_ot': r"^CT(?:\/GA)?\/CT:\s+(\d+)\s+\(\(converted\) top strand\)$",
                'aln_strand_ctot': r"^GA(?:\/CT)?\/CT:\s+(\d+)\s+\(complementary to \(converted\) top strand\)$",
                'aln_strand_ctob': r"^GA(?:\/CT)?\/GA:\s+(\d+)\s+\(complementary to \(converted\) bottom strand\)$",
                'aln_strand_ob': r"^CT(?:\/GA)?\/GA:\s+(\d+)\s+\(\(converted\) bottom strand\)$",
                'aln_strand_directional': r"^Option '--(directional)' specified \(default mode\): alignments to complementary strands \(CTOT, CTOB\) were ignored \(i.e. not performed\)$"
            },
            'dedup': {
                # 'aligned_reads' overwrites previous, but I trust this more
                # Leave the number from the alignment report in case deduplication is not run
                'aligned_reads': r"^Total number of alignments analysed in .+:\s+(\d+)$",
                'dup_reads': r"^Total number duplicated alignments removed:\s+(\d+)",
                'dup_reads_percent': r"^Total number duplicated alignments removed:\s+\d+\s+\(([\d\.]+)%\)",
                'dedup_reads': r"^Total count of deduplicated leftover sequences:\s+(\d+)",
                'dedup_reads_percent': r"^Total count of deduplicated leftover sequences:\s+\d+\s+\(([\d\.]+)% of total\)"
            },
            'methextract': {
                # These calls are typically done after deduplication
                'me_total_c': r"^Total number of C's analysed:\s+(\d+)$",
                'me_meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
                'me_meth_chg': r"^Total methylated C's in CHG context:\s+(\d+)",
                'me_meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
                'me_unmeth_cpg': r"^Total C to T conversions in CpG context:\s+(\d+)",
                'me_unmeth_chg': r"^Total C to T conversions in CHG context:\s+(\d+)",
                'me_unmeth_chh': r"^Total C to T conversions in CHH context:\s+(\d+)",
                'me_percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
                'me_percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
                'me_percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%"
            }
        }
        for sn, data in self.bismark_raw_data.items():
            for report_type in regexes.keys():
                for k, r in regexes[report_type].items():
                    try:
                        r_search = re.search(r, data[report_type], re.MULTILINE)
                        if r_search:
                            try:
                                self.bismark_data[sn][k] = float(r_search.group(1))
                            except ValueError:
                                self.bismark_data[sn][k] = r_search.group(1) # NaN
                    except KeyError:
                        pass # Missing report type

    def bismark_stats_table(self, report):
        """ Take the parsed stats from the Bismark reports and add them to the
        basic stats table at the top of the report """

        # Use several try blocks in case one of the report types is missing
        # If exception is triggered, header rows won't be added
        try:
            for sn, data in self.bismark_data.items():
                report['general_stats']['rows'][sn]['percent_cpg_meth'] = '<td class="text-right">{:.1f}%</td>'.format(data['me_percent_cpg_meth'])
                report['general_stats']['rows'][sn]['total_c'] = '<td class="text-right">{:.1f}</td>'.format(data['me_total_c']/1000000)
            report['general_stats']['headers']['percent_cpg_meth'] = '<th class="chroma-col" data-chroma-scale="BrBG" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: % Cytosines methylated in CpG context (meth&nbsp;extraction)">%&nbsp;Meth</span></th>'
            report['general_stats']['headers']['total_c'] = '<th class="chroma-col" data-chroma-scale="Purples" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total number of C\'s analysed, in millions (meth&nbsp;extraction)">M&nbsp;C\'s</span></th>'
        except KeyError:
            # Use numbers from alignment instead
            try:
                for sn, data in self.bismark_data.items():
                    report['general_stats']['rows'][sn]['percent_cpg_meth'] = '<td class="text-right">{:.1f}%</td>'.format(data['aln_percent_cpg_meth'])
                    report['general_stats']['rows'][sn]['total_c'] = '<td class="text-right">{:.1f}</td>'.format(data['aln_total_c']/1000000)
                report['general_stats']['headers']['percent_cpg_meth'] = '<th class="chroma-col" data-chroma-scale="Greens" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: % Cytosines methylated in CpG context (alignment)">%&nbsp;Meth</span></th>'
                report['general_stats']['headers']['total_c'] = '<th class="chroma-col" data-chroma-scale="Purples" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total number of C\'s analysed, in millions (alignment)">M&nbsp;C\'s</span></th>'
            except KeyError:
                pass

        try:
            for sn, data in self.bismark_data.items():
                report['general_stats']['rows'][sn]['bismark_dedup_reads_percent'] = '<td class="text-right">{:.1f}%</td>'.format(data['dup_reads_percent'])
                report['general_stats']['rows'][sn]['bismark_dedup_reads'] = '<td class="text-right">{:.1f}</td>'.format(data['dedup_reads']/1000000)
                report['general_stats']['rows'][sn]['bismark_aligned'] = '<td class="text-right">{:.1f}</td>'.format(data['aligned_reads']/1000000)
            report['general_stats']['headers']['bismark_dedup_reads_percent'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Percent Duplicated Alignments">%&nbsp;Dups</span></th>'
            report['general_stats']['headers']['bismark_dedup_reads'] = '<th class="chroma-col" data-chroma-scale="Greens" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Deduplicated Alignments (millions)">M&nbsp;Unique</span></th>'
        except KeyError:
            pass

        try:
            for sn, data in self.bismark_data.items():
                report['general_stats']['rows'][sn]['bismark_percent_aligned'] = '<td class="text-right">{:.1f}%</td>'.format((data['aligned_reads']/data['total_reads'])*100)
                report['general_stats']['rows'][sn]['bismark_aligned'] = '<td class="text-right">{:.1f}</td>'.format(data['aligned_reads']/1000000)
            report['general_stats']['headers']['bismark_percent_aligned'] = '<th class="chroma-col" data-chroma-scale="YlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Percent Aligned Sequences">%&nbsp;Aligned</span></th>'
            report['general_stats']['headers']['bismark_aligned'] = '<th class="chroma-col" data-chroma-scale="PuRd" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total Aligned Sequences (millions)">M&nbsp;Aligned</span></th>'
        except KeyError:
            pass


    def parse_alignment_chart_data (self):
        """ Make a data structure suitable for HighCharts for the alignment plot """
        self.bismark_sn_categories = list()
        series = OrderedDict()
        series['Deduplicated Unique Alignments'] = list()
        series['Duplicated Unique Alignments'] = list()
        series['Aligned Uniquely'] = list()
        series['Aligned Ambiguously'] = list()
        series['Did Not Align'] = list()
        series['No Genomic Sequence'] = list()
        colours = {
            'No Genomic Sequence': '#f28f43',
            'Did Not Align': '#0d233a',
            'Aligned Ambiguously': '#492970',
            'Aligned Uniquely': '#2f7ed8',
            'Duplicated Unique Alignments': '#2f7ed8',
            'Deduplicated Unique Alignments': '#8bbc21',
        }
        optional_cats = ['Aligned Uniquely', 'Duplicated Unique Alignments', 'Deduplicated Unique Alignments']
        for sn in sorted(self.bismark_data.keys()):
            self.bismark_sn_categories.append(sn)
            series['No Genomic Sequence'].append(int(self.bismark_data[sn].get('discarded_reads', 0)))
            series['Did Not Align'].append(int(self.bismark_data[sn].get('no_alignments', 0)))
            series['Aligned Ambiguously'].append(int(self.bismark_data[sn].get('ambig_reads', 0)))
            try:
                series['Duplicated Unique Alignments'].append(int(self.bismark_data[sn]['dup_reads']))
                series['Deduplicated Unique Alignments'].append(int(self.bismark_data[sn]['dedup_reads']))
                series['Aligned Uniquely'].append(0)
            except KeyError:
                series['Aligned Uniquely'].append(int(self.bismark_data[sn].get('aligned_reads', 0)))

        self.bismark_aln_plot_series = list()
        for cat in series:
            if cat not in optional_cats or (len(series[cat]) > 0 and max(series[cat]) > 0):
                self.bismark_aln_plot_series.append({
                    'name': cat,
                    'color': colours[cat],
                    'data': series[cat]
                })

    def bismark_alignment_chart (self):
        """ Make the HighCharts HTML to plot the alignment rates """

        return '<div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#bismark_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_percent" data-target="#bismark_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="bismark_alignment_plot" class="hc-plot"></div> \n\
        <script type="text/javascript"> \n\
            bismark_alignment_cats = {};\n\
            bismark_alignment_data = {};\n\
            var bismark_alignment_pconfig = {{ \n\
                "title": "Bismark Alignment Scores",\n\
                "ylab": "# Reads",\n\
                "ymin": 0,\n\
                "stacking": "normal" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bismark_alignment_plot", bismark_alignment_cats, bismark_alignment_data, bismark_alignment_pconfig); \
            }}); \
        </script>'.format(json.dumps(self.bismark_sn_categories), json.dumps(self.bismark_aln_plot_series));


    def parse_methylation_chart_data (self):
        """ Make a data structure suitable for HighCharts for the methylation plot """
        self.bismark_meth_helptext = "Numbers taken from methylation extraction report."
        self.bismark_meth_snames = list()
        series = {
            'cpg': OrderedDict([('Methylated CpG', list()), ('Unmethylated CpG', list())]),
            'chg': OrderedDict([('Methylated CHG', list()), ('Unmethylated CHG', list())]),
            'chh': OrderedDict([('Methylated CHH', list()), ('Unmethylated CHH', list())]),
        }
        colours = {
            'cpg': {'Methylated CpG': '#0d233a', 'Unmethylated CpG': '#2f7ed8'},
            'chg': {'Methylated CHG': '#1aadce', 'Unmethylated CHG': '#8bbc21'},
            'chh': {'Methylated CHH': '#492970', 'Unmethylated CHH': '#910000'}
        }
        for sn in sorted(self.bismark_data.keys()):
            self.bismark_meth_snames.append(sn)
            try:
                series['cpg']['Methylated CpG'].append(int(self.bismark_data[sn]['me_meth_cpg']))
                series['cpg']['Unmethylated CpG'].append(int(self.bismark_data[sn]['me_unmeth_cpg']))
                series['chg']['Methylated CHG'].append(int(self.bismark_data[sn]['me_meth_chg']))
                series['chg']['Unmethylated CHG'].append(int(self.bismark_data[sn]['me_unmeth_chg']))
                series['chh']['Methylated CHH'].append(int(self.bismark_data[sn]['me_meth_chh']))
                series['chh']['Unmethylated CHH'].append(int(self.bismark_data[sn]['me_unmeth_chh']))
            except KeyError:
                series['cpg']['Methylated CpG'].append(int(self.bismark_data[sn]['aln_meth_cpg']))
                series['cpg']['Unmethylated CpG'].append(int(self.bismark_data[sn]['aln_unmeth_cpg']))
                series['chg']['Methylated CHG'].append(int(self.bismark_data[sn]['aln_meth_chg']))
                series['chg']['Unmethylated CHG'].append(int(self.bismark_data[sn]['aln_unmeth_chg']))
                series['chh']['Methylated CHH'].append(int(self.bismark_data[sn]['aln_meth_chh']))
                series['chh']['Unmethylated CHH'].append(int(self.bismark_data[sn]['aln_unmeth_chh']))
                self.bismark_meth_helptext = "Numbers taken from Bismark alignment report"

        self.bismark_meth_cpg_series = list()
        self.bismark_meth_chg_series = list()
        self.bismark_meth_chh_series = list()
        for cat in series['cpg']:
            self.bismark_meth_cpg_series.append({
                'name': cat,
                'color': colours['cpg'][cat],
                'data': series['cpg'][cat]
            })
        for cat in series['chg']:
            self.bismark_meth_chg_series.append({
                'name': cat,
                'color': colours['chg'][cat],
                'data': series['chg'][cat]
            })
        for cat in series['chh']:
            self.bismark_meth_chh_series.append({
                'name': cat,
                'color': colours['chh'][cat],
                'data': series['chh'][cat]
            })

    def bismark_methlyation_chart (self):
        """ Make the HighCharts HTML to plot the methylation calls """

        return '<p class="text-muted">{}<p> \n\
        <div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm" data-action="set_numbers" data-target="#bismark_methylation_plot">Number of Calls</button> \n\
			<button class="btn btn-default btn-sm active" data-action="set_percent" data-target="#bismark_methylation_plot">Percentages</button> \n\
		</div> &nbsp; &nbsp; \n\
        <div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_data" data-newdata="bismark_methylation_cpg_data" data-target="#bismark_methylation_plot">CpG</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_data" data-newdata="bismark_methylation_chg_data" data-target="#bismark_methylation_plot">CHG</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_data" data-newdata="bismark_methylation_chh_data" data-target="#bismark_methylation_plot">CHH</button> \n\
		</div> \n\
        <div id="bismark_methylation_plot" class="hc-plot"></div> \n\
        <script type="text/javascript"> \n\
            bismark_methylation_cats = {};\n\
            bismark_methylation_cpg_data = {};\n\
            bismark_methylation_chg_data = {};\n\
            bismark_methylation_chh_data = {};\n\
            var bismark_methylation_pconfig = {{ \n\
                "title": "Cytosine Methylation",\n\
                "ylab": "% Calls",\n\
                "ymin": 0,\n\
                "stacking": "percent" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bismark_methylation_plot", bismark_methylation_cats, bismark_methylation_cpg_data, bismark_methylation_pconfig); \
            }}); \
        </script>'.format(self.bismark_meth_helptext, json.dumps(self.bismark_meth_snames), json.dumps(self.bismark_meth_cpg_series), json.dumps(self.bismark_meth_chg_series), json.dumps(self.bismark_meth_chh_series));


    def parse_strand_chart_data (self):
        """ Make a data structure suitable for HighCharts for the strand alignment plot """
        self.bismark_strand_samples = list()
        self.bismark_directional_mode = 0
        series = OrderedDict()
        series['OB'] = list()
        series['CTOB'] = list()
        series['CTOT'] = list()
        series['OT'] = list()
        for sn in sorted(self.bismark_data.keys()):
            self.bismark_strand_samples.append(sn)
            series['OB'].append(int(self.bismark_data[sn].get('aln_strand_ob', 0)))
            series['CTOB'].append(int(self.bismark_data[sn].get('aln_strand_ctob', 0)))
            series['CTOT'].append(int(self.bismark_data[sn].get('aln_strand_ctot', 0)))
            series['OT'].append(int(self.bismark_data[sn].get('aln_strand_ot', 0)))
            if 'aln_strand_directional' in self.bismark_data[sn]:
                self.bismark_directional_mode += 1

        if self.bismark_directional_mode == len(self.bismark_strand_samples):
            series.pop('CTOB', None)
            series.pop('CTOT', None)

        self.bismark_strand_plot_series = list()
        for cat in series:
            self.bismark_strand_plot_series.append({
                'name': cat,
                'data': series[cat]
            })

    def bismark_strand_chart (self):
        """ Make the HighCharts HTML to plot the strand alignment rates """
        if self.bismark_directional_mode == len(self.bismark_strand_samples):
            d_mode = '<p>All samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.</p>'
        elif self.bismark_directional_mode == 0:
            d_mode = ''
        else:
            d_mode = '<p>{} samples were run with <code>--directional</code> mode; alignments to complementary strands (CTOT, CTOB) were ignored.</p>'.format(self.bismark_directional_mode)

        return '{} <div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm" data-action="set_numbers" data-target="#bismark_strand_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm active" data-action="set_percent" data-target="#bismark_strand_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="bismark_strand_alignment_plot" class="hc-plot"></div>\n\
        <div class="row"><div class="col-sm-6"><ul class="list-unstyled">\n\
            <li><strong>OT</strong>: Original top strand</li>\n\
            <li><strong>OB</strong>: Original bottom strand</li>\n\
        </ul></div><div class="col-sm-6"><ul class="list-unstyled">\n\
            <li><strong>CTOT</strong>: Complementary to original top strand</li>\n\
            <li><strong>CTOB</strong>: Complementary to original bottom strand</li>\n\
        </ul></div></div>\n\
        <script type="text/javascript"> \n\
            bismark_strand_cats = {};\n\
            bismark_strand_data = {};\n\
            var bismark_strand_pconfig = {{ \n\
                "title": "Alignment to Individual Bisulfite Strands",\n\
                "ylab": "% Reads",\n\
                "ymin": 0,\n\
                "stacking": "percent" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bismark_strand_alignment_plot", bismark_strand_cats, bismark_strand_data, bismark_strand_pconfig); \
            }}); \
        </script>'.format(d_mode, json.dumps(self.bismark_strand_samples), json.dumps(self.bismark_strand_plot_series));
