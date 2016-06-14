#!/usr/bin/env python

""" MultiQC Submodule to parse output from Qualimap BamQC """

from __future__ import print_function
import logging
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Qualimap BamQC reports and parse their data """
    
    sp = config.sp['qualimap']['bamqc']
    
    # General stats - genome_results.txt
    self.qualimap_bamqc_genome_results = dict()
    for f in self.find_log_files(sp['genome_results']):
        parse_genome_results(self, f)
    
    # Coverage - coverage_histogram.txt
    self.qualimap_bamqc_coverage_hist = dict()
    for f in self.find_log_files(sp['coverage'], filehandles=True):
        parse_coverage(self, f)
    
    # Insert size - insert_size_histogram.txt
    self.qualimap_bamqc_insert_size_hist = dict()
    for f in self.find_log_files(sp['insert_size'], filehandles=True):
        parse_insert_size(self, f)
    
    # Genome fraction - genome_fraction_coverage.txt
    self.qualimap_bamqc_genome_fraction_cov = dict()
    for f in self.find_log_files(sp['genome_fraction'], filehandles=True):
        parse_genome_fraction(self, f)
    
    # GC distribution - mapped_reads_gc-content_distribution.txt
    self.qualimap_bamqc_gc_content_dist = dict()
    for f in self.find_log_files(sp['gc_dist'], filehandles=True):
        parse_gc_dist(self, f)
    
    # Make the plots for the report
    report_sections(self)
    
    # Set up the general stats table
    general_stats_headers(self)
    
    # Return the number of reports we found
    return len(self.qualimap_bamqc_genome_results.keys())
    
def parse_genome_results(self, f):
    """ Parse the contents of the Qualimap BamQC genome_results.txt file """
    regexes = {
        'bam_file': r"bam file = (.+)",
        'total_reads': r"number of reads = ([\d,]+)",
        'mapped_reads': r"number of mapped reads = ([\d,]+)",
        'mapped_bases': r"number of mapped bases = ([\d,]+)",
        'sequenced_bases': r"number of sequenced bases = ([\d,]+)",
        'mean_insert_size': r"mean insert size = ([\d,\.]+)",
        'median_insert_size': r"median insert size = ([\d,\.]+)",
        'mean_mapping_quality': r"mean mapping quality = ([\d,\.]+)",
    }
    d = dict()
    for k, r in regexes.items():
        r_search = re.search(r, f['f'], re.MULTILINE)
        if r_search:
            try:
                d[k] = float(r_search.group(1).replace(',',''))
            except ValueError:
                d[k] = r_search.group(1)
    
    # Check we have an input filename
    if 'bam_file' not in d:
        log.debug("Couldn't find an input filename in genome_results file {}".format(f['fn']))
        return None
    
    # Get a nice sample name
    s_name = self.clean_s_name(d['bam_file'], f['root'])
    
    # Add to general stats table & calculate a nice % aligned
    try:
        self.general_stats_data[s_name]['total_reads'] = d['total_reads']
        self.general_stats_data[s_name]['mapped_reads'] = d['mapped_reads']
        d['percentage_aligned'] = (d['mapped_reads'] / d['total_reads'])*100
        self.general_stats_data[s_name]['percentage_aligned'] = d['percentage_aligned']
    except KeyError:
        pass
    
    # Save results
    if s_name in self.qualimap_bamqc_genome_results:
        log.debug("Duplicate genome results sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_genome_results[s_name] = d
    self.add_data_source(f, s_name=s_name, section='genome_results')
    
    
def parse_coverage(self, f):
    """ Parse the contents of the Qualimap BamQC Coverage Histogram file """
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/coverage_histogram.txt
    s_name = self.get_s_name(f)
    
    d = dict()
    for l in f['f']:
        if l.startswith('#'):
            continue
        coverage, count = l.split(None, 1)
        coverage = int(round(float(coverage)))
        count = float(count)
        d[coverage] = count
    
    if len(d) == 0:
        log.debug("Couldn't parse contents of coverage histogram file {}".format(f['fn']))
        return None
    
    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    median_coverage = None
    for thiscov, thiscount in d.items():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_coverage = thiscov
            break
    self.general_stats_data[s_name]['median_coverage'] = median_coverage
    
    # Save results
    if s_name in self.qualimap_bamqc_coverage_hist:
        log.debug("Duplicate coverage histogram sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_coverage_hist[s_name] = d
    self.add_data_source(f, s_name=s_name, section='coverage_histogram')
                
def parse_insert_size(self, f):
    """ Parse the contents of the Qualimap BamQC Insert Size Histogram file """
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/insert_size_histogram.txt
    s_name = self.get_s_name(f)
    
    d = dict()
    zero_insertsize = 0
    for l in f['f']:
        if l.startswith('#'):
            continue
        insertsize, count = l.split(None, 1)
        insertsize = int(round(float(insertsize)))
        count = float(count) / 1000000
        if(insertsize == 0):
            zero_insertsize = count
        else:
            d[insertsize] = count
    
    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    median_insert_size = None
    for thisins, thiscount in d.items():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_insert_size = thisins
            break
    # Add the median insert size to the general stats table
    self.general_stats_data[s_name]['median_insert_size'] = median_insert_size
    
    # Save results
    if s_name in self.qualimap_bamqc_insert_size_hist:
        log.debug("Duplicate insert size histogram sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_insert_size_hist[s_name] = d
    self.add_data_source(f, s_name=s_name, section='insert_size_histogram')
                
def parse_genome_fraction(self, f):
    """ Parse the contents of the Qualimap BamQC Genome Fraction Coverage file """
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/genome_fraction_coverage.txt
    s_name = self.get_s_name(f)
    
    d = dict()
    fifty_x_pc = thirty_x_pc = ten_x_pc = five_x_pc = one_x_pc = 100
    for l in f['f']:
        if l.startswith('#'):
            continue
        try:
            coverage, percentage = l.split(None, 1)
        except ValueError:
            continue
        coverage = int(round(float(coverage)))
        percentage = float(percentage)
        d[coverage] = percentage
        if coverage <= 50 and fifty_x_pc > percentage:
            fifty_x_pc = percentage
        if coverage <= 30 and thirty_x_pc > percentage:
            thirty_x_pc = percentage
        if coverage <= 10 and ten_x_pc > percentage:
            ten_x_pc = percentage
        if coverage <= 5 and five_x_pc > percentage:
            five_x_pc = percentage
        if coverage <= 1 and one_x_pc > percentage:
            one_x_pc = percentage
    
    # Add the coverage cutoffs to the general stats table
    self.general_stats_data[s_name]['fifty_x_pc'] = fifty_x_pc
    self.general_stats_data[s_name]['thirty_x_pc'] = thirty_x_pc
    self.general_stats_data[s_name]['ten_x_pc'] = ten_x_pc
    self.general_stats_data[s_name]['five_x_pc'] = five_x_pc
    self.general_stats_data[s_name]['one_x_pc'] = one_x_pc
    
    # Save results
    if s_name in self.qualimap_bamqc_genome_fraction_cov:
        log.debug("Duplicate genome fraction coverage sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_genome_fraction_cov[s_name] = d
    self.add_data_source(f, s_name=s_name, section='genome_fraction_coverage')
                
def parse_gc_dist(self, f):
    """ Parse the contents of the Qualimap BamQC Mapped Reads GC content distribution file """
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt
    s_name = self.get_s_name(f)
    
    d = dict()
    avg_gc = 0
    for l in f['f']:
        if l.startswith('#'):
            continue
        sections = l.split(None, 2)
        gc = int(round(float(sections[0])))
        content = float(sections[1])
        avg_gc += gc * content
        d[gc] = content
    
    # Add average GC to the general stats table
    self.general_stats_data[s_name]['avg_gc'] = avg_gc
    
    # Save results
    if s_name in self.qualimap_bamqc_gc_content_dist:
        log.debug("Duplicate Mapped Reads GC content distribution sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_gc_content_dist[s_name] = d
    self.add_data_source(f, s_name=s_name, section='mapped_gc_distribution')


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

def general_stats_headers (self):
    self.general_stats_headers['avg_gc'] = {
        'title': 'Avg. GC',
        'description': 'Average GC content',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'Set1',
        'format': '{:.0f}%'
    }
    self.general_stats_headers['median_insert_size'] = {
        'title': 'Insert Size',
        'description': 'Median Insert Size',
        'min': 0,
        'suffix': 'bp',
        'scale': 'PuOr',
        'format': '{:.0f}'
    }
    self.general_stats_headers['fifty_x_pc'] = {
        'title': '&ge; 50X',
        'description': 'Fraction of genome with at least 50X coverage',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn',
        'format': '{:.1f}%',
        'hidden': True
    }
    self.general_stats_headers['thirty_x_pc'] = {
        'title': '&ge; 30X',
        'description': 'Fraction of genome with at least 30X coverage',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn',
        'format': '{:.1f}%'
    }
    self.general_stats_headers['ten_x_pc'] = {
        'title': '&ge; 10X',
        'description': 'Fraction of genome with at least 10X coverage',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn',
        'format': '{:.1f}%',
        'hidden': True
    }
    self.general_stats_headers['five_x_pc'] = {
        'title': '&ge; 05X',
        'description': 'Fraction of genome with at least 05X coverage',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn',
        'format': '{:.1f}%',
        'hidden': True
    }
    self.general_stats_headers['one_x_pc'] = {
        'title': '&ge; 01X',
        'description': 'Fraction of genome with at least 01X coverage',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'RdYlGn',
        'format': '{:.1f}%',
        'hidden': True
    }
    self.general_stats_headers['median_coverage'] = {
        'title': 'Coverage',
        'description': 'Median coverage',
        'min': 0,
        'suffix': 'X',
        'scale': 'BuPu'
    }
    self.general_stats_headers['percentage_aligned'] = {
        'title': '% Aligned',
        'description': '% mapped reads',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'YlGn',
        'format': '{:.1f}%'
    }
    self.general_stats_headers['mapped_reads'] = {
        'title': 'Aligned',
        'description': 'Number of mapped reads (millions)',
        'min': 0,
        'scale': 'RdYlGn',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
        'hidden': True
    }
    self.general_stats_headers['total_reads'] = {
        'title': 'Total Reads',
        'description': 'Number of reads (millions)',
        'min': 0,
        'scale': 'Blues',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
        'hidden': True
    }
