#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Submodule to handle code for Qualimap BamQC """

from __future__ import print_function
from collections import OrderedDict
import io
import logging
import os

from collections import defaultdict

from multiqc import config, BaseMultiqcModule

def parse_reports(self):
    """ Find Qualimap BamQC reports and parse their data """
    
    self.qualimap_bamqc_coverage_hist = dict()
    self.qualimap_bamqc_insert_size_hist = dict()
    self.qualimap_bamqc_genome_fraction_cov = dict()
    self.qualimap_bamqc_gc_content_dist = dict()
    
    # Find QualiMap reports
    qualimap_raw_data = {}
    for directory in config.analysis_dir:
        for root, dirnames, filenames in os.walk(directory, followlinks=True):
            raw_data_dir = 'raw_data'
            for d in dirnames:
                if raw_data_dir in d:
                    raw_data_dir = d
            if 'genome_results.txt' in filenames and raw_data_dir in dirnames:
                with io.open(os.path.join(root, 'genome_results.txt'), 'r') as gr:
                    for l in gr:
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root)
                
                if s_name in qualimap_raw_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        
                #### Coverage histogram
                cov_report = os.path.join(root, raw_data_dir, 'coverage_histogram.txt')
                if os.path.exists(cov_report):
                    self.qualimap_bamqc_coverage_hist[s_name] = {}
                    with io.open(cov_report, 'r') as fh:
                        next(fh)
                        for l in fh:
                            coverage, count = l.split(None, 1)
                            coverage = int(round(float(coverage)))
                            count = float(count)
                            self.qualimap_bamqc_coverage_hist[s_name][coverage] = count
                    # Find median
                    num_counts = sum(self.qualimap_bamqc_coverage_hist[s_name].values())
                    cum_counts = 0
                    median_coverage = None
                    for thiscov, thiscount in self.qualimap_bamqc_coverage_hist[s_name].items():
                        cum_counts += thiscount
                        if cum_counts >= num_counts/2:
                            median_coverage = thiscov
                            break
                    # Add median to the general stats table
                    self.general_stats[s_name]['median_coverage'] = median_coverage
                
                
                #### Insert size histogram
                ins_size = os.path.join(root, raw_data_dir, 'insert_size_histogram.txt')
                if os.path.exists(ins_size):
                    self.qualimap_bamqc_insert_size_hist[s_name] = {}
                    zero_insertsize = 0
                    with io.open(ins_size, 'r') as fh:
                        next(fh)
                        for l in fh:
                            insertsize, count = l.split(None, 1)
                            insertsize = int(round(float(insertsize)))
                            count = float(count) / 1000000
                            if(insertsize == 0):
                                zero_insertsize = count
                            else:
                                self.qualimap_bamqc_insert_size_hist[s_name][insertsize] = count
                    # Find median
                    num_counts = sum(self.qualimap_bamqc_insert_size_hist[s_name].values())
                    cum_counts = 0
                    median_insert_size = None
                    for thisins, thiscount in self.qualimap_bamqc_insert_size_hist[s_name].items():
                        cum_counts += thiscount
                        if cum_counts >= num_counts/2:
                            median_insert_size = thisins
                            break
                    # Add the median insert size to the general stats table
                    self.general_stats[s_name]['median_insert_size'] = median_insert_size
                
                
                #### Genome Fraction Coverage
                frac_cov = os.path.join(root, raw_data_dir, 'genome_fraction_coverage.txt')
                if os.path.exists(frac_cov):
                    self.qualimap_bamqc_genome_fraction_cov[s_name] = {}
                    thirty_x_pc = 100
                    max_obs_x = 0
                    halfway_cov = None
                    with io.open(frac_cov, 'r') as fh:
                        next(fh)
                        for l in fh:
                            coverage, percentage = l.split(None, 1)
                            coverage = int(round(float(coverage)))
                            percentage = float(percentage)
                            self.qualimap_bamqc_genome_fraction_cov[s_name][coverage] = percentage
                            if coverage <= 30 and thirty_x_pc > percentage:
                                thirty_x_pc = percentage
                    # Add the median % genome >= 30X coverage to the general stats table
                    self.general_stats[s_name]['thirty_x_pc'] = thirty_x_pc
                
                
                #### GC Distribution
                gc_report = os.path.join(root, raw_data_dir, 'mapped_reads_gc-content_distribution.txt')
                if os.path.exists(gc_report):
                    self.qualimap_bamqc_gc_content_dist[s_name] = {}
                    avg_gc = 0
                    with io.open(gc_report, 'r') as fh:
                        next(fh)
                        for l in fh:
                            sections = l.split(None, 2)
                            gc = int(round(float(sections[0])))
                            content = float(sections[1])
                            avg_gc += gc * content
                            self.qualimap_bamqc_gc_content_dist[s_name][gc] = content
                    # Add average GC to the general stats table
                    self.general_stats[s_name]['avg_gc'] = avg_gc


def report_sections(self):
    """ Add results from Qualimap BamQC parsing to the report """
    # Append to self.sections list
    
    # Section 1 - BamQC Coverage Histogram
    if len(self.qualimap_bamqc_coverage_hist) > 0:
        # Chew back on histogram to prevent long flat tail
        # (find a sensible max x - lose 1% of longest tail)
        max_x = 0
        for d in self.qualimap_bamqc_coverage_hist.values():
            total = sum(d.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / total > 0.01:
                    max_x = max(max_x, count)
                    break                    
        self.sections.append({
            'name': 'Coverage Histogram',
            'anchor': 'qualimap-coverage-histogram',
            'content': self.plot_xy_data(self.qualimap_bamqc_coverage_hist, {
                'title': 'Coverage Histogram',
                'ylab': 'Genome Bin Counts',
                'xlab': 'Coverage (X)',
                'ymin': 0,
                'xmin': 0,
                'xmax': max_x,
                'xDecimals': False,
                'tt_label': '<b>{point.x}X</b>: {point.y}',
            })
        })

    # Section 2 - Insert size histogram
    if len(self.qualimap_bamqc_insert_size_hist) > 0:
        self.sections.append({
            'name': 'Insert size Histogram',
            'anchor': 'qualimap-insert-size-histogram',
            'content': self.plot_xy_data(self.qualimap_bamqc_insert_size_hist, {
                'title': 'Insert Size Histogram',
                'ylab': 'Fraction of reads',
                'xlab': 'Insert Size (bp)',
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            })
        })

    # Section 3 - Genome Fraction coverage
    if len(self.qualimap_bamqc_genome_fraction_cov) > 0:
        self.sections.append({
            'name': 'Genome Fraction Coverage',
            'anchor': 'qualimap-genome-fraction-coverage',
            'content': self.plot_xy_data(self.qualimap_bamqc_genome_fraction_cov, {
                'title': 'Genome Fraction Coverage',
                'ylab': 'Fraction of reference (%)',
                'xlab': 'Coverage (X)',
                'ymax': 100,
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
            })
        })

    # Section 4 - GC-content distribution
    if len(self.qualimap_bamqc_gc_content_dist) > 0:
        self.sections.append({
            'name': 'GC-content distribution',
            'anchor': 'qualimap-gc-distribution',
            'content': self.plot_xy_data(self.qualimap_bamqc_gc_content_dist, {
                'title': 'GC-content distribution',
                'ylab': 'Fraction of reads',
                'xlab': 'GC content (%)',
                'ymin': 0,
                'xmin': 0,
                'xmax': 100,
                'tt_label': '<b>{point.x}%</b>: {point.y:.3f}',
            })
        })



def stats_table(self):
    """ Take the parsed stats from the QualiMap report and add them to the
    basic stats table at the top of the report """
    
    headers = OrderedDict()
    headers['median_coverage'] = {
        'title': 'Coverage',
        'description': 'Median coverage',
        'min': 0,
        'scale': 'RdBu'
    }
    headers['median_insert_size'] = {
        'title': 'Insert Size',
        'description': 'Median Insert Size',
        'min': 0,
        'scale': 'PuOr',
        'format': '{:.0f}'
    }
    headers['thirty_x_pc'] = {
        'title': '&ge; 30X',
        'description': 'Fraction of genome with at least 30X coverage',
        'max': 100,
        'min': 0,
        'scale': 'RdYlGn',
        'format': '{:.1f}%'
    }
    headers['avg_gc'] = {
        'title': 'Avg. GC',
        'description': 'Average GC content',
        'max': 80,
        'min': 20,
        'format': '{:.0f}%'
    }
    self.general_stats_addcols(self.general_stats, headers)

