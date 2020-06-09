#!/usr/bin/env python
from __future__ import print_function

import copy
import re
import os
import json
from collections import OrderedDict, defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.dragen.utils import Metric
from multiqc.plots import linegraph, bargraph, heatmap, table, boxplot
from multiqc.utils import report

# Initialise the logger
import logging

log = logging.getLogger(__name__)

N_QV = 2
ADAPTER_SEQS = ['AGATCGGAAGAG', 'ATGGAATTCTCG', 'CTGTCTCTTATA']


class DragenFastQcMetrics(BaseMultiqcModule):
    """
    Rendering ALL THE THINGS!
    """

    def add_fastqc_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/fastqc_metrics'):
            data = parse_fastqc_metrics_file(f)
            self.add_data_source(f, section='stats')
            data_by_sample.update(data)

        # Filter to strip out ignored sample names:
        self.fastqc_data = self.ignore_samples(data_by_sample)
        if not self.fastqc_data:
            return
        log.info('Found time metrics for {} samples'.format(len(self.fastqc_data)))

        # Now add each section in order
        self.positional_quality_range_plot()
        self.positional_mean_quality_plot()
        self.per_seq_quality_plot()
        self.n_content_plot()
        self.gc_content_plot()
        self.gc_content_mean_quality_plot()
        self.seq_length_dist_plot()
        self.sequence_content_plot()
        self.adapter_content_plot()

        return self.fastqc_data.keys()

    def positional_quality_range_plot(self):
        """STUFF"""

        data = OrderedDict()
        GROUP = "POSITIONAL QUALITY"

        for s_name in sorted(self.fastqc_data):
            data[s_name] = defaultdict(float)

            sorted_keys = sortPosQualTableKeys(self.fastqc_data[s_name][GROUP])
            for key in sorted_keys:
                value = int(self.fastqc_data[s_name][GROUP][key])
                parts = key.split()
                pos = average_from_range(parts[1])
                quantile = int(parts[2][:-1])
                qv = int(value)

                try:
                    data[s_name][pos][quantile] = qv
                except:
                    data[s_name][pos] = OrderedDict()
                    data[s_name][pos][quantile] = qv

        pconfig = {
            'id': 'fastqc_per_base_sequence_quality_plot',
            'title': 'DRAGEN-QC: Per-Position Quality Range',
            'ylab': 'Phred Quality Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'ymax': 43,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            # 'colors': self.get_status_cols('per_base_sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }

        self.add_section(
            name='Per-Position Quality Score Ranges',
            anchor='fastqc_pos_qual_ranges',
            description='The range of quality value across each base position in each sample or read',
            plot=boxplot.plot(data, pconfig)
        )

    def positional_mean_quality_plot(self):
        """ Create the HTML for the positional mean-quality score plot """

        base_dict = {"A": {}, "C": {}, "G": {}, "T": {}, "N": {}}

        AVG_GROUP = "POSITIONAL BASE MEAN QUALITY"
        COUNT_GROUP = "POSITIONAL BASE CONTENT"

        data = dict()
        for s_name in self.fastqc_data:
            # Parse our per-base, per-position average qualities into a dictionary
            avgs = copy.deepcopy(base_dict)
            for key, value in self.fastqc_data[s_name][AVG_GROUP].items():
                if value == "NA":
                    continue
                parts = key.split()
                pos = average_from_range(parts[1])
                base = parts[2].upper()
                avgs[base][pos] = float(value)

            # Parse matching per-base and total counts by position
            counts = copy.deepcopy(base_dict)
            totals = defaultdict(int)
            for key, value in self.fastqc_data[s_name][COUNT_GROUP].items():
                parts = key.split()
                pos = average_from_range(parts[1])
                base = parts[2].upper()
                counts[base][pos] = int(value)
                totals[pos] += int(value)

            # Use the the count and averages to recompute total QVs
            qv_sums = defaultdict(int)
            for base, pos_data in counts.items():
                for pos, count in pos_data.items():
                    if count == 0:
                        continue
                    elif base == "N":
                        qv_sums[pos] += count * N_QV
                    else:
                        qv_sums[pos] += int(round(count * avgs[base][pos]))

            # Compute the positional, base-agnostic mean QV
            data[s_name] = dict()
            for pos, qv_sum in qv_sums.items():
                total = totals[pos]
                if total > 0:
                    data[s_name][int(pos)] = qv_sum / total

        pconfig = {
            'id': 'fastqc_per_base_sequence_quality_plot',
            'title': 'DRAGEN-QC: Per-Position Quality Scores',
            'ylab': 'Phred Quality Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            # 'colors': self.get_status_cols('per_base_sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }

        self.add_section(
            name='Per-Position Mean Quality Scores',
            anchor='fastqc_per_base_sequence_quality',
            description='The mean quality value across each base position in the read.',
            helptext='''
            To enable multiple samples to be plotted on the same graph, only the mean quality
            scores are plotted (unlike the box plots seen in FastQC reports).

            Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):

            _The y-axis on the graph shows the quality scores. The higher the score, the better
            the base call. The background of the graph divides the y axis into very good quality
            calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
            The quality of calls on most platforms will degrade as the run progresses, so it is
            common to see base calls falling into the orange area towards the end of a read._
            ''',
            plot=linegraph.plot(data, pconfig)
        )

    def per_seq_quality_plot(self):
        """ Create the HTML for the per sequence quality score plot """

        data = dict()
        GROUP = "READ MEAN QUALITY"
        MAX_QV = 64
        max_non_zero = 0
        for s_name in self.fastqc_data:
            data[s_name] = dict()
            group_data = self.fastqc_data[s_name][GROUP]
            for qv in range(MAX_QV):
                metric = "Q{0} Reads".format(qv)
                count = group_data[metric]
                if count > 0:
                    max_non_zero = max(qv, max_non_zero)
                data[s_name][qv] = count

        for s_name in data:
            for qv in range(max_non_zero+2, MAX_QV):
                del data[s_name][qv]

        if len(data) == 0:
            log.debug('per_seq_quality not found in DRAGEN FastQC reports')
            return None

        pconfig = {
            'id': 'dragenqc_per_sequence_quality_scores_plot',
            'title': 'DRAGEN-QC: Per-Sequence Quality Scores',
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Quality Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            # 'colors': self.get_status_cols('per_sequence_quality_scores'),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.add_section(
            name='Per-Sequence Quality Scores',
            anchor='dragenqc_per_sequence_quality_scores',
            description='The number of reads with average quality scores. Shows if a subset of reads has poor quality.',
            helptext='''
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html):
            _The per sequence quality score report allows you to see if a subset of your
            sequences have universally low quality values. It is often the case that a
            subset of sequences will have universally poor quality, however these should
            represent only a small percentage of the total sequences._
            ''',
            plot=linegraph.plot(data, pconfig)
        )

    def n_content_plot(self):
        """ Create the HTML for the per base N content plot """

        data = dict()
        totals = defaultdict(int)
        non_n = defaultdict(int)
        GROUP = "POSITIONAL BASE CONTENT"
        for s_name in self.fastqc_data:
            # Count total bases
            total_group_data = self.fastqc_data[s_name][GROUP]
            for metric, value in total_group_data.items():
                avg_pos = average_pos_from_metric(metric)
                totals[avg_pos] += value
                base = metric.split()[2]
                if base != "N":
                    non_n[avg_pos] += value

            # Convert Total and Non-N counts into N%
            data[s_name] = dict()
            for pos, count in totals.items():
                if count == 0:
                    continue
                non_n_count = non_n[pos]
                n_count = count - non_n_count
                n_frac = 100.0 * n_count / float(count)
                data[s_name][pos] = n_frac

        if len(data) == 0:
            log.debug('per_base_n_content not found in DRAGEN FastQC reports')
            return None

        pconfig = {
            'id': 'dragenqc_per_base_n_content_plot',
            'title': 'DRAGEN-QC: Per-Position N Content',
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            # 'colors': self.get_status_cols('per_base_n_content'),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }

        self.add_section(
            name='Per-Position N Content',
            anchor='dragenqc_per_base_n_content',
            description='The percentage of base calls at each position for which an `N` was called.',
            helptext='''
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html):
            _If a sequencer is unable to make a base call with sufficient confidence then it will
            normally substitute an `N` rather than a conventional base call. This graph shows the
            percentage of base calls at each position for which an `N` was called._
            _It's not unusual to see a very low proportion of Ns appearing in a sequence, especially
            nearer the end of a sequence. However, if this proportion rises above a few percent
            it suggests that the analysis pipeline was unable to interpret the data well enough to
            make valid base calls._
            ''',
            plot=linegraph.plot(data, pconfig)
        )

    def gc_content_plot(self):
        """ Create the HTML for the FastQC GC content plot """

        data = dict()
        data_norm = dict()
        LEN_GROUP = "READ LENGTHS"
        GC_GROUP = "READ GC CONTENT"
        for s_name in self.fastqc_data:
            data[s_name] = defaultdict(float)

            # First figure out the baseline
            max_len = 0
            for metric, value in self.fastqc_data[s_name][LEN_GROUP].items():
                if int(value) > 0:
                    pos = average_from_range(metric.split()[0][:-2])
                    max_len = max(pos, max_len)

            group_data = self.fastqc_data[s_name][GC_GROUP]
            for metric, value in self.fastqc_data[s_name][GC_GROUP].items():
                pct = percentage_from_content_metric(metric)
                data[s_name][pct] = value

            data_norm[s_name] = dict()
            total = sum([c for c in data[s_name].values()])
            for gc, count in data[s_name].items():
                if total > 0:
                    data_norm[s_name][gc] = (count / total) * 100

        if len(data) == 0:
            log.debug('per_sequence_gc_content not found in FastQC reports')
            return None

        pconfig = {
            'id': 'dragenqc_per_sequence_gc_content_plot',
            'title': 'DRAGEN-QC: Per-Sequence GC Content',
            'xlab': '% GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            # 'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'Count'}
            ]
        }

        self.add_section(
            name='Per-Sequence GC Content',
            anchor='dragenqc_per_sequence_gc_content',
            description='''The average GC content of reads. Normal random library typically have a
            roughly normal distribution of GC content.''',
            helptext='''
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html):
            _This module measures the GC content across the whole length of each sequence
            in a file and compares it to a modelled normal distribution of GC content._
            _In a normal random library you would expect to see a roughly normal distribution
            of GC content where the central peak corresponds to the overall GC content of
            the underlying genome. Since we don't know the the GC content of the genome the
            modal GC content is calculated from the observed data and used to build a
            reference distribution._
            _An unusually shaped distribution could indicate a contaminated library or
            some other kinds of biased subset. A normal distribution which is shifted
            indicates some systematic bias which is independent of base position. If there
            is a systematic bias which creates a shifted normal distribution then this won't
            be flagged as an error by the module since it doesn't know what your genome's
            GC content should be._
            ''',
            plot=linegraph.plot([data_norm, data], pconfig)
        )

    def gc_content_mean_quality_plot(self):
        """ Create the HTML for the positional mean-quality score plot """

        GROUP = "READ GC CONTENT QUALITY"

        data = dict()
        for s_name in self.fastqc_data:
            data[s_name] = dict()

            for key, value in self.fastqc_data[s_name][GROUP].items():
                parts = key.split()
                pct = int(parts[0][:-1])
                try:
                    data[s_name][pct] = float(value)
                except:
                    continue

        pconfig = {
            'id': 'fastqc_gc_content_mean_sequence_quality_plot',
            'title': 'DRAGEN-QC: GC Content Mean Quality Scores',
            'ylab': 'Phred Quality Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            # 'colors': self.get_status_cols('per_base_sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }

        self.add_section(
            name='GC Content Mean Quality Scores',
            anchor='fastqc_gc_content_mean_sequence_quality',
            description='The mean quality value across each base position in the read.',
            helptext='''
            To enable multiple samples to be plotted on the same graph, only the mean quality
            scores are plotted (unlike the box plots seen in FastQC reports).

            Taken from the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html):

            _The y-axis on the graph shows the quality scores. The higher the score, the better
            the base call. The background of the graph divides the y axis into very good quality
            calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
            The quality of calls on most platforms will degrade as the run progresses, so it is
            common to see base calls falling into the orange area towards the end of a read._
            ''',
            plot=linegraph.plot(data, pconfig)
        )

    def seq_length_dist_plot(self):
        """ Create the HTML for the Sequence Length Distribution plot """

        data = dict()
        seq_lengths = set()
        multiple_lenths = False
        avg_to_range = dict()
        GROUP = "READ LENGTHS"
        for s_name in self.fastqc_data:
            data[s_name] = dict()

            group_data = self.fastqc_data[s_name][GROUP]
            for metric, value in group_data.items():
                if value > 0:
                    avg_pos = average_pos_from_size(metric)
                    data[s_name][avg_pos] = value
                    avg_to_range[avg_pos] = metric.split('bp')[0]

            seq_lengths.update([avg_to_range[k] for k in data[s_name].keys()])

            if len(set(data[s_name].keys())) > 1:
                multiple_lenths = True

        if len(data) == 0:
            log.debug('sequence_length_distribution not found in FastQC reports')
            return None

        if not multiple_lenths:
            lengths = 'bp , '.join([str(l) for l in list(seq_lengths)])
            desc = 'All samples have sequences within a single length bin ({}bp).'.format(lengths)
            if len(seq_lengths) > 1:
                desc += ' See the <a href="#general_stats">General Statistics Table</a>.'
            self.add_section(
                name='Sequence Length Distribution',
                anchor='dragenqc_sequence_length_distribution',
                description='<div class="alert alert-info">{}</div>'.format(desc)
            )
        else:
            pconfig = {
                'id': 'dragenqc_sequence_length_distribution_plot',
                'title': 'DRAGEN-QC: Sequence Length Distribution',
                'ylab': 'Read Count',
                'xlab': 'Sequence Length (bp)',
                'ymin': 0,
                'yMinTickInterval': 0.1,
                'xDecimals': False,
                # 'colors': self.get_status_cols('sequence_length_distribution'),
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            }
            self.add_section(
                name='Sequence Length Distribution',
                anchor='dragenqc_sequence_length_distribution',
                description='''The distribution of fragment sizes (read lengths) found.
                    See the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)''',
                plot=linegraph.plot(data, pconfig)
            )

    def sequence_content_plot(self):
        """ Create the epic HTML for the FastQC sequence content heatmap """

        # Prep the data
        data = dict()
        GROUP = "POSITIONAL BASE CONTENT"
        for s_name in sorted(self.fastqc_data.keys()):
            data[s_name] = dict()
            group_data = self.fastqc_data[s_name][GROUP]

            totals = defaultdict(int)
            for metric, value in group_data.items():
                parts = metric.split()
                #avg_pos = average_from_range(parts[1])
                avg_pos = parts[1]
                base = parts[-2].lower()

                if avg_pos not in data[s_name]:
                    data[s_name][avg_pos] = dict()

                # Store the current count and add it to the total
                data[s_name][avg_pos][base] = value
                totals[avg_pos] += value

            # Use the accumulated totals to normalize each bin to a percentage
            for pos, total in totals.items():
                if total == 0:
                    del data[s_name][pos]
                    continue
                for base in "acgt":
                    try:
                        data[s_name][pos][base] = (float(data[s_name][pos][base])/float(total)) * 100.0
                    except:
                        pass
                data[s_name][pos]["base"] = pos

        if len(data) == 0:
            log.debug('sequence_content not found in FastQC reports')
            return None

        html = '''<div id="fastqc_per_base_sequence_content_plot_div">
            <div class="alert alert-info">
               <span class="glyphicon glyphicon-hand-up"></span>
               Click a sample row to see a line plot for that dataset.
            </div>
            <h5><span class="s_name text-primary"><span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name</span></h5>
            <button id="fastqc_per_base_sequence_content_export_btn"><span class="glyphicon glyphicon-download-alt"></span> Export Plot</button>
            <div class="fastqc_seq_heatmap_key">
                Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
                <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
                <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
            </div>
            <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
                <div id="{id}" class="fastqc_per_base_sequence_content_plot hc-plot has-custom-export">
                    <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
                </div>
            </div>
            <div class="clearfix"></div>
        </div>
        <script type="application/json" class="fastqc_seq_content">{d}</script>
        '''.format(
            # Generate unique plot ID, needed in mqc_export_selectplots
            id=report.save_htmlid('fastqc_per_base_sequence_content_plot'),
            d=json.dumps([self.anchor.replace('-', '_'), data]),
        )

        self.add_section(
            name='Per-Position Sequence Content',
            anchor='fastqc_per_base_sequence_content',
            description='The proportion of each base position for which each of the four normal DNA bases has been called.',
            helptext='''
            To enable multiple samples to be shown in a single plot, the base composition data
            is shown as a heatmap. The colours represent the balance between the four bases:
            an even distribution should give an even muddy brown colour. Hover over the plot
            to see the percentage of the four bases under the cursor.
            **To see the data as a line plot, as in the original FastQC graph, click on a sample track.**
            From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html):
            _Per Base Sequence Content plots out the proportion of each base position in a
            file for which each of the four normal DNA bases has been called._
            _In a random library you would expect that there would be little to no difference
            between the different bases of a sequence run, so the lines in this plot should
            run parallel with each other. The relative amount of each base should reflect
            the overall amount of these bases in your genome, but in any case they should
            not be hugely imbalanced from each other._
            _It's worth noting that some types of library will always produce biased sequence
            composition, normally at the start of the read. Libraries produced by priming
            using random hexamers (including nearly all RNA-Seq libraries) and those which
            were fragmented using transposases inherit an intrinsic bias in the positions
            at which reads start. This bias does not concern an absolute sequence, but instead
            provides enrichement of a number of different K-mers at the 5' end of the reads.
            Whilst this is a true technical bias, it isn't something which can be corrected
            by trimming and in most cases doesn't seem to adversely affect the downstream
            analysis._
            ''',
            content=html
        )

    def adapter_content_plot(self):
        """ Create the epic HTML for the FastQC adapter content plot"""

        # Prep the data
        data = dict()
        COUNT_GROUP = "POSITIONAL BASE CONTENT"
        ADP_GROUP = "SEQUENCE POSITIONS"
        for s_name in sorted(self.fastqc_data.keys()):

            totals = defaultdict(int)
            for key, value in self.fastqc_data[s_name][COUNT_GROUP].items():
                parts = key.split()
                pos = average_from_range(parts[1])
                totals[pos] += int(value)

            adps = defaultdict(int)
            for key, value in self.fastqc_data[s_name][ADP_GROUP].items():
                parts = key.split()
                seq = parts[0].split("'")[1]
                if seq not in ADAPTER_SEQS:
                    continue
                pos = average_from_range(parts[1][:-2])
                adps[pos] += int(value)

            data[s_name] = dict()
            cumsum = 0
            for pos, adp_count in sorted(adps.items()):
                total = totals[pos]
                cumsum += adp_count
                if total > 0 and cumsum > 0:
                    data[s_name][pos] = 100.0 * cumsum / total

        pconfig = {
            'id': 'fastqc_adapter_content_plot',
            'title': 'FastQC: Adapter Content',
            'ylab': '% of Sequences',
            'xlab': 'Position (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'hide_empty': True,
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ],
        }

        self.add_section(
            name='Adapter Content',
            anchor='fastqc_adapter_content',
            description='''The cumulative percentage count of the proportion of your
            library which has seen each of the adapter sequences at each position.''',
            helptext='''
            Note that only samples with &ge; 0.1% adapter contamination are shown.

            There may be several lines per sample, as one is shown for each adapter
            detected in the file.

            From the [FastQC Help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html):

            _The plot shows a cumulative percentage count of the proportion
            of your library which has seen each of the adapter sequences at each position.
            Once a sequence has been seen in a read it is counted as being present
            right through to the end of the read so the percentages you see will only
            increase as the read length goes on._
            ''',
            plot=linegraph.plot(data, pconfig)
        )


def parse_fastqc_metrics_file(f):
    """
    NA12878.fastqc_metrics.csv

    READ MEAN QUALITY,Read1,Q30 Reads,151104
    ...
    POSITIONAL BASE CONTENT,Read1,ReadPos 1 A Bases,2963930
    ...
    NUCLEOTIDE QUALITY,Read1,Q7 A Bases,777933
    ...
    READ LENGTHS,Read1,145-152bp Length Reads,10142922
    ...
    READ BASE CONTENT,Read1,0% A Reads,1997
    ...
    """
    f['s_name'] = re.search(r'(.*).fastqc_metrics.csv', f['fn']).group(1)
    r1_name = "{}_R1".format(f['s_name'])
    r2_name = "{}_R2".format(f['s_name'])

    data = {}
    data[r1_name] = defaultdict(lambda: defaultdict(int))
    data[r2_name] = defaultdict(lambda: defaultdict(int))
    for line in f['f'].splitlines():
        group, mate, metric, value = line.split(',')
        try:
            value = int(value)
        except ValueError:
            pass

        # Store each value by group and by metric
        if mate == "Read1":
            s_name = r1_name
        elif mate == "Read2":
            s_name = r2_name

        data[s_name][group][metric] = value

    # Delete empty mate groups so we don't generate empty datasets
    for s_name in [r1_name, r2_name]:
        if len(data[s_name]) == 0:
            del data[s_name]

    return data


def average_from_range(metric_range):
    if metric_range.startswith('>='):
        metric_range = metric_range[2:]
    if metric_range.endswith('+'):
        metric_range = metric_range[:-1]
    if "-" in metric_range:
        start, end = metric_range.split('-')
        avg_pos = (int(end) + int(start)) / 2.0
    else:
        avg_pos = int(metric_range)
    return avg_pos


def average_pos_from_metric(metric):
    parts = metric.split()
    metric_range = parts[1]
    return average_from_range(metric_range)


def average_pos_from_size(metric):
    parts = metric.split()
    len_range = parts[0].split('bp')[0]
    return average_from_range(len_range)


def percentage_from_content_metric(metric):
    parts = metric.split()
    pct = int(parts[0].split("%")[0])
    return pct


def base_from_content_metric(metric):
    parts = metric.split()
    return parts[1]


def pos_qual_table_cmp(key):
    parts = key.split()
    pos = average_from_range(parts[1])
    pct = int(parts[2][:-1])

    return (pos * 1000 + pct)


def sortPosQualTableKeys(data_dict):
    return sorted(data_dict.keys(), key=pos_qual_table_cmp)


def generateGcModel(read_length):
    models = {}

    claimed = [0] * 101
    for i in range(read_length):
        pct = round(100 * i / read_length)
        claimed[pct+1] += 1

    return claimed
