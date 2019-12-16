#!/usr/bin/env python
from __future__ import print_function

import itertools
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
import logging
log = logging.getLogger(__name__)


class DragenCoverage(BaseMultiqcModule):
    def parse_coverage(self):
        all_data_by_sample = dict()

        for f in itertools.chain(self.find_log_files('dragen/wgs_contig_mean_cov_normal'),
                                 self.find_log_files('dragen/wgs_contig_mean_cov_tumor')):
            data_by_sample = parse_wgs_contig_mean_cov(f)
            if data_by_sample:
                for sn, data in data_by_sample.items():
                    if sn in all_data_by_sample:
                        log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
                    self.add_data_source(f, section='stats')

                    all_data_by_sample[sn] = data

        # Filter to strip out ignored sample names:
        all_data_by_sample = self.ignore_samples(all_data_by_sample)

        if not all_data_by_sample:
            return
        log.info('Found variant calling metrics for {} samples'.format(len(all_data_by_sample)))

        # self.add_section(
        #     name='Fragment length histogram',
        #     anchor='dragen-fragment-length-histogram',
        #     description='Distribution of estimated fragment lengths of mapped reads.',
        #     plot=linegraph.plot(all_data_by_sample, {
        #         'id': 'dragen_fragment_length',
        #         'title': 'Dragen: fragment length histogram',
        #         'ylab': 'Fraction of reads',
        #         'xlab': 'Fragment length (bp)',
        #         'ymin': 0,
        #         'xmin': 0,
        #         'tt_label': '<b>{point.x} bp</b>: {point.y}',
        #     })
        # )


def parse_wgs_contig_mean_cov(f):
    """
    T_SRR7890936_50pc.wgs_contig_mean_cov_normal.csv
    T_SRR7890936_50pc.wgs_contig_mean_cov_tumor.csv

    chr1,11292297134,48.9945
    chr10,6482885699,48.6473
    ...
    chrUn_GL000218v1,20750824,128.77
    chrX,3590295769,23.1792
    chrY,42229820,1.5987
    chrY_KI270740v1_random,0,0
    Autosomal regions ,130912665915,47.4953
    """

    data_by_sample = defaultdict(dict)

    # sample = None
    # for line in f['f'].splitlines():
    #     if line.startswith('#Sample'):
    #         sample = line.split('#Sample: ')[1]
    #     else:
    #         assert sample is not None
    #         frag_len, cnt = line.split(',')
    #         try:
    #             frag_len = int(frag_len)
    #             cnt = int(cnt)
    #         except ValueError:
    #             assert line == 'FragmentLength,Count', line
    #         else:
    #             data_by_sample[sample][frag_len] = cnt

    return data_by_sample


