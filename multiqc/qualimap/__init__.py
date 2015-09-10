#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
import io
import json
import logging
import os

from collections import defaultdict, OrderedDict

import multiqc

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Qualimap'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "QualiMap"
        self.anchor = "qualimap"
        self.intro = '<p><a href="http://qualimap.bioinfo.cipf.es/" target="_blank">QualiMap</a> \
             is a platform-independent application to facilitate the quality control of alignment \
             sequencing data and its derivatives like feature counts.</p>'
        self.analysis_dir = report['analysis_dir']
        self.parsed_stats = defaultdict(dict)

        # Find QualiMap reports
        qualimap_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            if 'genome_results.txt' in filenames and 'raw_data_qualimapReport' in dirnames:
                with io.open(os.path.join(root, 'genome_results.txt'), 'r') as gr:
                    for l in gr:
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]), root, prepend_dirs=report['prepend_dirs'])

                s_name = self.clean_s_name(s_name, root, prepend_dirs=report['prepend_dirs'])
                if s_name in qualimap_raw_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

                qualimap_raw_data[s_name] = {}
                qualimap_raw_data[s_name]['reports'] = {os.path.splitext(r)[0]: os.path.join(root, 'raw_data_qualimapReport', r) \
                    for r in os.listdir(os.path.join(root, 'raw_data_qualimapReport'))}

        if len(qualimap_raw_data) == 0:
            log.debug("Could not find any reports in {}".format(self.analysis_dir))
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

        # General stats table
        self.qualimap_stats_table(report)


    def qualimap_cov_his(self, qualimap_raw_data):
        parsed_data = {}
        for sn, data in qualimap_raw_data.iteritems():
            cov_report = data['reports'].get('coverage_histogram')
            if cov_report:
                counts=OrderedDict()
                try:
                    with io.open(cov_report, 'r') as fh:
                        next(fh) # skip the header
                        for line in fh:
                            (coverage, count) = line.split(None, 1)
                            coverage = int(round(float(coverage)))
                            count = float(count)
                            counts[coverage] = count
                except IOError as e:
                    log.error("Could not load input file: {}".format(cov_report))
                    raise

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
                counts = OrderedDict()
                zero_insertsize = 0
                try:
                    with open(ins_size, 'r') as fh:
                        next(fh) # skip the header
                        for line in fh:
                            (insertsize, count) = line.split(None, 1)
                            insertsize = int(round(float(insertsize)))
                            count = float(count) / 1000000
                            if(insertsize == 0):
                                zero_insertsize = count
                            else:
                                counts[insertsize] = count
                except IOError as e:
                    logging.error("Could not load input file: {}".format(fn))
                    raise

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


    def qualimap_stats_table(self, report):
        """ Take the parsed stats from the QualiMap report and add them to the
        basic stats table at the top of the report """

        # General stats table headers
        report['general_stats']['headers']['median_coverage'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Qualimap: Median coverage">Med. Cov</span></th>'
        report['general_stats']['headers']['median_insert_size'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Qualimap: Median Insert Size">Med. Ins</span></th>'

        rowcounts = { 'median_coverage' : 0, 'median_insert_size': 0}

        for samp, vals in self.parsed_stats.items():
            for k, v in vals.items():
                report['general_stats']['rows'][samp][k] = '<td class="text-right">{}</td>'.format(v)
                rowcounts[k] += 1

        # Remove header if we don't have any filled cells for it
        for k in rowcounts.keys():
            if rowcounts[k] == 0:
                report['general_stats']['headers'].pop(k, None)
