#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
import io
import json
import logging
import os

from collections import defaultdict

import multiqc
from multiqc import config

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Qualimap'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(log)

        # Static variables
        self.name = "QualiMap"
        self.anchor = "qualimap"
        self.intro = '<p><a href="http://qualimap.bioinfo.cipf.es/" target="_blank">QualiMap</a> \
             is a platform-independent application to facilitate the quality control of alignment \
             sequencing data and its derivatives like feature counts.</p>'

        self.parsed_stats = defaultdict(dict)

        # Find QualiMap reports
        qualimap_raw_data = {}
        for root, dirnames, filenames in os.walk(config.analysis_dir, followlinks=True):
            raw_data_dir = 'raw_data'
            for d in dirnames:
                if raw_data_dir in d:
                    raw_data_dir = d
            if 'genome_results.txt' in filenames and raw_data_dir in dirnames:
                with io.open(os.path.join(root, 'genome_results.txt'), 'r') as gr:
                    for l in gr:
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root)

                s_name = self.clean_s_name(s_name, root)
                if s_name in qualimap_raw_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

                qualimap_raw_data[s_name] = {}
                qualimap_raw_data[s_name]['reports'] = {os.path.splitext(r)[0]: os.path.join(root, raw_data_dir, r) \
                    for r in os.listdir(os.path.join(root, raw_data_dir))}

        if len(qualimap_raw_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(qualimap_raw_data)))

        self.sections = list()

        # Section 1 - Coverage Histogram
        histogram_data = self.qualimap_cov_his(qualimap_raw_data)
        if len(histogram_data) > 0:
            self.sections.append({
                'name': 'Coverage Histogram',
                'anchor': 'qualimap-coverage-histogram',
                'content': self.plot_xy_data(histogram_data, {
                    'title': 'Coverage Histogram',
                    'ylab': 'Genome Bin Counts',
                    'xlab': 'Coverage (X)',
                    'ymin': 0,
                    'xmin': 0,
                    'tt_label': '<b>{point.x}-X coverage </b>',
                })
            })

        # Section 2 - Insert size histogram
        histogram_data = self.qualimap_ins_size_his(qualimap_raw_data)
        if len(histogram_data) > 0:
            self.sections.append({
                'name': 'Insert size Histogram',
                'anchor': 'qualimap-insert-size-histogram',
                'content': self.plot_xy_data(histogram_data, {
                    'title': 'Insert Size Histogram',
                    'ylab': 'Number of reads',
                    'xlab': 'Insert Size (bp)',
                    'ymin': 0,
                    'xmin': 0,
                    'tt_label': '<b>{point.x} insert size (bd) </b>',
                })
            })

        # Section 3 - Genome Fraction coverage
        histogram_data = self.qualimap_gen_frac_his(qualimap_raw_data)
        if len(histogram_data) > 0:
            self.sections.append({
                'name': 'Genome Fraction Coverage',
                'anchor': 'qualimap-genome-fraction-coverage',
                'content': self.plot_xy_data(histogram_data, {
                    'title': 'Genome Fraction Coverage',
                    'ylab': 'Fraction of reference (%)',
                    'xlab': 'Coverage (X)',
                    'ymin': 0,
                    'xmin': 0,
                    'tt_label': '<b>{point.x}-X coverage </b>',
                })
            })

            # Section 4 - GC-content distribution
            histogram_data = self.qualimap_gc_distribution(qualimap_raw_data)
            if len(histogram_data) > 0:
                self.sections.append({
                    'name': 'GC-content distribution',
                    'anchor': 'qualimap-gc-distribution',
                    'content': self.plot_xy_data(histogram_data, {
                        'title': 'GC-content distribution',
                        'ylab': 'Fraction of reads',
                        'xlab': 'GC content (%)',
                        'ymin': 0,
                        'xmin': 0,
                        'tt_label': '<b>GC-content (%) </b>',
                    })
                })

        # General stats table
        self.qualimap_stats_table()


    def qualimap_gc_distribution(self, qualimap_raw_data):
        parsed_data = {}
        for sn, data in qualimap_raw_data.iteritems():
            gc_report = data['reports']['mapped_reads_gc-content_distribution']
            if gc_report:
                counts={}
                avg_gc = 0
                with io.open(gc_report, 'r') as fh:
                    next(fh)
                    for l in fh:
                        sections = l.split(None, 2)
                        gc = int(round(float(sections[0])))
                        cont = float(sections[1])
                        avg_gc += gc*cont
                        counts[gc] = cont

                parsed_data[sn] = counts

                #Add reads avg. GC to the general stats table
                self.parsed_stats[sn]['avg_gc'] = avg_gc

        return parsed_data

    def qualimap_cov_his(self, qualimap_raw_data):
        parsed_data = {}
        for sn, data in qualimap_raw_data.iteritems():
            cov_report = data['reports'].get('coverage_histogram')
            if cov_report:
                counts={}
                with io.open(cov_report, 'r') as fh:
                    next(fh)
                    for l in fh:
                        coverage, count = l.split(None, 1)
                        coverage = int(round(float(coverage)))
                        count = float(count)
                        counts[coverage] = count

                parsed_data[sn] = counts

                # Find median
                num_counts = sum(counts.values())
                cum_counts = 0
                median_coverage = None
                for thiscov, thiscount in counts.iteritems():
                    cum_counts += thiscount
                    if cum_counts >= num_counts/2:
                        median_coverage = thiscov
                        break

                # Add median to the general stats table
                self.parsed_stats[sn]['median_coverage'] = median_coverage

        return parsed_data


    def qualimap_ins_size_his(self, qualimap_raw_data):
        parsed_data = {}
        for sn, data in qualimap_raw_data.iteritems():
            ins_size = data['reports'].get('insert_size_histogram')
            if ins_size:
                counts = {}
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
                            counts[insertsize] = count

                parsed_data[sn] = counts

                # Find median
                num_counts = sum(counts.values())
                cum_counts = 0
                median_insert_size = None
                for thisins, thiscount in counts.iteritems():
                    cum_counts += thiscount
                    if cum_counts >= num_counts/2:
                        median_insert_size = thisins
                        break

                # Add the median insert size to the general stats table
                self.parsed_stats[sn]['median_insert_size'] = median_insert_size

        return parsed_data


    def qualimap_gen_frac_his(self, qualimap_raw_data):
        parsed_data = {}
        for sn, data in qualimap_raw_data.iteritems():
            frac_cov = data['reports'].get('genome_fraction_coverage')
            if frac_cov:
                thirty_x_pc = 100
                max_obs_x = 0
                halfway_cov = None
                counts={}
                with io.open(frac_cov, 'r') as fh:
                    next(fh)
                    for l in fh:
                        coverage, percentage = l.split(None, 1)
                        coverage = int(round(float(coverage)))
                        percentage = float(percentage)
                        counts[coverage] = percentage

                        if coverage <= 30 and thirty_x_pc > percentage:
                            thirty_x_pc = percentage

                parsed_data[sn] = counts

                # Add the median % genome >= 30X coverage to the general stats table
                self.parsed_stats[sn]['thirty_x_pc'] = thirty_x_pc

        return parsed_data


    def qualimap_stats_table(self):
        """ Take the parsed stats from the QualiMap report and add them to the
        basic stats table at the top of the report """

        # General stats table headers
        config.general_stats['headers']['median_coverage'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev"><span data-toggle="tooltip" title="Qualimap: Median coverage">Coverage</span></th>'
        config.general_stats['headers']['median_insert_size'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev"><span data-toggle="tooltip" title="Qualimap: Median Insert Size">Insert Size</span></th>'
        config.general_stats['headers']['thirty_x_pc'] = '<th class="chroma-col" data-chroma-scale="RdYlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Qualimap: Genome Fraction Coverage">% G Covered</span></th>'
        config.general_stats['headers']['avg_gc'] = '<th class="chroma-col" data-chroma-scale="BrBG" data-chroma-max="80" data-chroma-min="20"><span data-toggle="tooltip" title="Qualimap: Average GC content">Avg. GC</span></th>'

        rowcounts = { 'median_coverage' : 0, 'median_insert_size': 0,
                      'thirty_x_pc': 0, 'avg_gc': 0}

        for samp, vals in self.parsed_stats.items():
            for k, v in vals.items():
                v = round(v, 3) if type(v) == float else v
                config.general_stats['rows'][samp][k] = '<td class="text-right">{}</td>'.format(v)
                rowcounts[k] += 1

        # Remove header if we don't have any filled cells for it
        for k in rowcounts.keys():
            if rowcounts[k] == 0:
                config.general_stats['headers'].pop(k, None)
