#!/usr/bin/env python

""" MultiQC Submodule to parse output from Qualimap BamQC """

from __future__ import print_function
import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Qualimap BamQC reports and parse their data """

    try:
        covs = config.qualimap_config['general_stats_coverage']
        assert type(covs) == list
        assert len(covs) > 0
        covs = [str(i) for i in covs]
        log.debug("Custom Qualimap thresholds: {}".format(", ".join([i for i in covs])))
    except (AttributeError, TypeError, AssertionError):
        covs = [1, 5, 10, 30, 50]
        covs = [str(i) for i in covs]
        log.debug("Using default Qualimap thresholds: {}".format(", ".join([i for i in covs])))
    self.covs = covs

    # General stats - genome_results.txt
    self.qualimap_bamqc_genome_results = dict()
    for f in self.find_log_files('qualimap/bamqc/genome_results'):
        parse_genome_results(self, f)
    self.qualimap_bamqc_genome_results = self.ignore_samples(self.qualimap_bamqc_genome_results)

    # Coverage - coverage_histogram.txt
    self.qualimap_bamqc_coverage_hist = dict()
    for f in self.find_log_files('qualimap/bamqc/coverage', filehandles=True):
        parse_coverage(self, f)
    self.qualimap_bamqc_coverage_hist = self.ignore_samples(self.qualimap_bamqc_coverage_hist)

    # Insert size - insert_size_histogram.txt
    self.qualimap_bamqc_insert_size_hist = dict()
    for f in self.find_log_files('qualimap/bamqc/insert_size', filehandles=True):
        parse_insert_size(self, f)
    self.qualimap_bamqc_insert_size_hist = self.ignore_samples(self.qualimap_bamqc_insert_size_hist)

    # GC distribution - mapped_reads_gc-content_distribution.txt
    self.qualimap_bamqc_gc_content_dist = dict()
    self.qualimap_bamqc_gc_by_species = dict()  # {'HUMAN': data_dict, 'MOUSE': data_dict}
    for f in self.find_log_files('qualimap/bamqc/gc_dist', filehandles=True):
        parse_gc_dist(self, f)
    self.qualimap_bamqc_gc_by_species = self.ignore_samples(self.qualimap_bamqc_gc_by_species)

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

def parse_gc_dist(self, f):
    """ Parse the contents of the Qualimap BamQC Mapped Reads GC content distribution file """
    # Get the sample name from the parent parent directory
    # Typical path: <sample name>/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt
    s_name = self.get_s_name(f)

    d = dict()
    reference_species = None
    reference_d = dict()
    avg_gc = 0
    for l in f['f']:
        if l.startswith('#'):
            sections = l.strip("\n").split("\t", 3)
            if len(sections) > 2:
                reference_species = sections[2]
            continue
        sections = l.strip("\n").split("\t", 3)
        gc = int(round(float(sections[0])))
        content = float(sections[1])
        avg_gc += gc * content
        d[gc] = content
        if len(sections) > 2:
            reference_content = float(sections[2])
            reference_d[gc] = reference_content

    # Add average GC to the general stats table
    self.general_stats_data[s_name]['avg_gc'] = avg_gc

    # Save results
    if s_name in self.qualimap_bamqc_gc_content_dist:
        log.debug("Duplicate Mapped Reads GC content distribution sample name found! Overwriting: {}".format(s_name))
    self.qualimap_bamqc_gc_content_dist[s_name] = d
    if reference_species and reference_species not in self.qualimap_bamqc_gc_by_species:
        self.qualimap_bamqc_gc_by_species[reference_species] = reference_d
    self.add_data_source(f, s_name=s_name, section='mapped_gc_distribution')


def report_sections(self):
    """ Add results from Qualimap BamQC parsing to the report """
    # Append to self.sections list

    if len(self.qualimap_bamqc_coverage_hist) > 0:
        # Chew back on histogram to prevent long flat tail
        # (find a sensible max x - lose 1% of longest tail)
        max_x = 0
        total_bases_by_sample = dict()
        for s_name, d in self.qualimap_bamqc_coverage_hist.items():
            total_bases_by_sample[s_name] = sum(d.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / total_bases_by_sample[s_name] > 0.01:
                    max_x = max(max_x, count)
                    break

        rates_within_threshs = dict()
        for s_name, hist in self.qualimap_bamqc_coverage_hist.items():
            total = total_bases_by_sample[s_name]
            rates_within_threshs[s_name] = _calculate_bases_within_thresholds(hist, total, range(max_x + 1))
            for c in self.covs:
                if int(c) in rates_within_threshs[s_name]:
                    self.general_stats_data[s_name]['{}_x_pc'.format(c)] = rates_within_threshs[s_name][int(c)]
                else:
                    self.general_stats_data[s_name]['{}_x_pc'.format(c)] = 0

        # Section 1 - BamQC Coverage Histogram
        self.add_section (
            name = 'Coverage histogram',
            anchor = 'qualimap-coverage-histogram',
            plot = linegraph.plot(self.qualimap_bamqc_coverage_hist, {
                'id': 'qualimap_coverage_histogram',
                'title': 'Coverage histogram',
                'ylab': 'Genome bin counts',
                'xlab': 'Coverage (X)',
                'ymin': 0,
                'xmin': 0,
                'xmax': max_x,
                'xDecimals': False,
                'tt_label': '<b>{point.x}X</b>: {point.y}',
            })
        )
        # Section 2 - BamQC cumulative coverage genome fraction
        self.add_section (
            name = 'Cumulative coverage genome fraction',
            anchor = 'qualimap-cumulative-genome-fraction-coverage',
            plot = linegraph.plot(rates_within_threshs, {
                'id': 'qualimap_genome_fraction',
                'title': 'Genome fraction covered by at least X reads',
                'ylab': 'Fraction of reference (%)',
                'xlab': 'Coverage (X)',
                'ymax': 100,
                'ymin': 0,
                'xmin': 0,
                'xmax': max_x,
                'xDecimals': False,
                'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
            })
        )

    # Section 3 - Insert size histogram
    if len(self.qualimap_bamqc_insert_size_hist) > 0:
        self.add_section (
            name = 'Insert size histogram',
            anchor = 'qualimap-insert-size-histogram',
            plot = linegraph.plot(self.qualimap_bamqc_insert_size_hist, {
                'id': 'qualimap_insert_size',
                'title': 'Insert size histogram',
                'ylab': 'Fraction of reads',
                'xlab': 'Insert Size (bp)',
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            })
        )

    # Section 4 - GC-content distribution
    if len(self.qualimap_bamqc_gc_content_dist) > 0:
        extra_series = []
        for i, (species_name, species_data) in enumerate(sorted(self.qualimap_bamqc_gc_by_species.items())):
            extra_series.append({
                'name': species_name,
                'data': list(species_data.items()),
                'dashStyle': 'Dash',
                'lineWidth': 1,
                'color': ['#000000', '#E89191'][i % 2],
            })
        desc = ''
        lg_config = {
            'id': 'qualimap_gc_content',
            'title': 'GC content distribution',
            'ylab': 'Fraction of reads',
            'xlab': 'GC content (%)',
            'ymin': 0,
            'xmin': 0,
            'xmax': 100,
            'tt_label': '<b>{point.x}%</b>: {point.y:.3f}'
        }
        if len(extra_series) == 1:
            desc += 'The dotted line represents a pre-calculated GC destribution for the reference genome.'
            lg_config['extra_series'] = extra_series
        elif len(extra_series) > 1:
            desc += 'The dotted lines represent pre-calculated GC destributions for the reference genomes.'
            lg_config['extra_series'] = extra_series

        self.add_section (
            name = 'GC content distribution',
            anchor = 'qualimap-gc-distribution',
            description = desc,
            plot = linegraph.plot(self.qualimap_bamqc_gc_content_dist, lg_config)
        )

def general_stats_headers (self):
    try:
        hidecovs = config.qualimap_config['general_stats_coverage_hidden']
        assert type(hidecovs) == list
        log.debug("Hiding Qualimap thresholds: {}".format(", ".join([i for i in hidecovs])))
    except (AttributeError, TypeError, AssertionError):
        hidecovs = [1, 5, 10, 50]
    hidecovs = [str(i) for i in hidecovs]

    self.general_stats_headers['avg_gc'] = {
        'title': 'Avg. GC',
        'description': 'Average GC content',
        'max': 100,
        'min': 0,
        'suffix': '%',
        'scale': 'Set1',
        'format': '{:,.0f}'
    }
    self.general_stats_headers['median_insert_size'] = {
        'title': 'Insert Size',
        'description': 'Median insert size',
        'min': 0,
        'suffix': 'bp',
        'scale': 'PuOr',
        'format': '{:,.0f}'
    }
    for c in self.covs:
        self.general_stats_headers['{}_x_pc'.format(c)] = {
            'title': '&ge; {}X'.format(c),
            'description': 'Fraction of genome with at least {}X coverage'.format(c),
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
            'hidden': c in hidecovs
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
        'scale': 'YlGn'
    }
    self.general_stats_headers['mapped_reads'] = {
        'title': '{} Aligned'.format(config.read_count_prefix),
        'description': 'Number of mapped reads ({})'.format(config.read_count_desc),
        'min': 0,
        'scale': 'RdYlGn',
        'shared_key': 'read_count',
        'modify': lambda x: x * config.read_count_multiplier,
        'hidden': True
    }
    self.general_stats_headers['total_reads'] = {
        'title': '{} Total reads'.format(config.read_count_prefix),
        'description': 'Number of reads ({})'.format(config.read_count_desc),
        'min': 0,
        'scale': 'Blues',
        'shared_key': 'read_count',
        'modify': lambda x: x * config.read_count_multiplier,
        'hidden': True
    }


def _calculate_bases_within_thresholds(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
    rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)

    for depth, bases in bases_by_depth.items():
        for t in depth_thresholds:
            if depth >= t:
                bases_within_threshs[t] += bases
    for t in depth_thresholds:
        bs = bases_within_threshs[t]
        if total_size > 0:
            rate = 100.0 * bases_within_threshs[t] / total_size
            assert rate <= 100, 'Error: rate is > 1: rate = ' + str(rate) + ', bases = ' + str(bs) + ', size = ' + str(total_size)
            rates_within_threshs[t] = rate
    return rates_within_threshs
