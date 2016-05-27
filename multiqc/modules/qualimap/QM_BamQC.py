#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Submodule to handle code for Qualimap BamQC """

import io
import logging
import os

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Qualimap BamQC reports and parse their data """
    
    self.qualimap_bamqc_coverage_hist = dict()
    self.qualimap_bamqc_insert_size_hist = dict()
    self.qualimap_bamqc_genome_fraction_cov = dict()
    self.qualimap_bamqc_gc_content_dist = dict()
    
    sp = config.sp['qualimap']['bamqc']
    
    # Find QualiMap reports
    for directory in config.analysis_dir:
        for root, dirnames, filenames in os.walk(directory, followlinks=True):
            raw_data_dir = sp['raw_data']
            for d in dirnames:
                if raw_data_dir in d:
                    raw_data_dir = d
            if sp['genome_results'] in filenames and raw_data_dir in dirnames:
                with io.open(os.path.join(root, sp['genome_results']), 'r') as gr:
                    for l in gr:
                        
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root)
                        
                        rhs = l.split(' = ')[-1].replace(',','')
                        num = rhs.split('(')[0]
                        if 'number of reads' in l:
                            self.general_stats[s_name]['total_reads'] = float( num )
                        if 'number of mapped reads' in l:
                            self.general_stats[s_name]['mapped_reads'] = float( num )
                
                # Calculate percentage aligned
                for s_name, v in self.general_stats.items():
                    if 'mapped_reads' in v and 'total_reads' in v:
                        self.general_stats[s_name]['percentage_aligned'] = (v['mapped_reads'] / v['total_reads'])*100
                
                if s_name in self.qualimap_bamqc_coverage_hist:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(s_name=s_name, source=os.path.abspath(os.path.join(root, sp['genome_results'])), section='genome_results')
        
                #### Coverage histogram
                cov_report = os.path.join(root, raw_data_dir, sp['coverage'])
                if os.path.exists(cov_report):
                    self.add_data_source(s_name=s_name, source=os.path.abspath(cov_report), section='coverage_histogram')
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
                ins_size = os.path.join(root, raw_data_dir, sp['insert_size'])
                if os.path.exists(ins_size):
                    self.add_data_source(s_name=s_name, source=os.path.abspath(ins_size), section='insert_size_histogram')
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
                frac_cov = os.path.join(root, raw_data_dir, sp['genome_fraction'])
                if os.path.exists(frac_cov):
                    self.add_data_source(s_name=s_name, source=os.path.abspath(frac_cov), section='genome_fraction')
                    self.qualimap_bamqc_genome_fraction_cov[s_name] = {}
                    fifty_x_pc = 100
                    thirty_x_pc = 100
                    ten_x_pc = 100
                    five_x_pc = 100
                    one_x_pc = 100
                    max_obs_x = 0
                    halfway_cov = None
                    with io.open(frac_cov, 'r') as fh:
                        next(fh)
                        for l in fh:
                            coverage, percentage = l.split(None, 1)
                            coverage = int(round(float(coverage)))
                            percentage = float(percentage)
                            self.qualimap_bamqc_genome_fraction_cov[s_name][coverage] = percentage
                            if coverage <= 50 and five_x_pc > percentage:
                                fifty_x_pc = percentage
                            if coverage <= 30 and thirty_x_pc > percentage:
                                thirty_x_pc = percentage
                            if coverage <= 10 and ten_x_pc > percentage:
                                ten_x_pc = percentage
                            if coverage <= 5 and five_x_pc > percentage:
                                five_x_pc = percentage
                            if coverage <= 1 and one_x_pc > percentage:
                                one_x_pc = percentage
                    # Add the median % genome >= 30X coverage to the general stats table
                    self.general_stats[s_name]['fifty_x_pc'] = fifty_x_pc
                    self.general_stats[s_name]['thirty_x_pc'] = thirty_x_pc
                    self.general_stats[s_name]['ten_x_pc'] = ten_x_pc
                    self.general_stats[s_name]['five_x_pc'] = five_x_pc
                    self.general_stats[s_name]['one_x_pc'] = one_x_pc
                
                
                #### GC Distribution
                gc_report = os.path.join(root, raw_data_dir, sp['gc_dist'])
                if os.path.exists(gc_report):
                    self.add_data_source(s_name=s_name, source=os.path.abspath(gc_report), section='gc_distribution')
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
            'content': plots.linegraph.plot(self.qualimap_bamqc_coverage_hist, {
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
            'content': plots.linegraph.plot(self.qualimap_bamqc_insert_size_hist, {
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
            'content': plots.linegraph.plot(self.qualimap_bamqc_genome_fraction_cov, {
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
            'content': plots.linegraph.plot(self.qualimap_bamqc_gc_content_dist, {
                'title': 'GC-content distribution',
                'ylab': 'Fraction of reads',
                'xlab': 'GC content (%)',
                'ymin': 0,
                'xmin': 0,
                'xmax': 100,
                'tt_label': '<b>{point.x}%</b>: {point.y:.3f}',
            })
        })

