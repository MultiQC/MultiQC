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

from .util import average_from_range, percentage_from_content_metric

# Initialise the logger
import logging

log = logging.getLogger(__name__)

N_QV = 2
ADAPTER_SEQS = ['AGATCGGAAGAG', 'ATGGAATTCTCG', 'CTGTCTCTTATA']


class DragenGcMetrics(BaseMultiqcModule):
    """
    Rendering ALL THE THINGS!
    """

    def add_gc_metrics(self):
        # Now add each section in order
        self.gc_content_plot()
        self.gc_content_mean_quality_plot()

        return self.fastqc_data.keys()

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
            'ylab': 'Count',
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
