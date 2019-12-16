#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenVCMetrics(BaseMultiqcModule):
    def parse_vc_metrics(self):
        all_data_by_sample = dict()

        for f in self.find_log_files('dragen/vc_metrics'):
            data_by_sample = parse_vc_metrics_file(f)
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


def parse_vc_metrics_file(f):
    """
    T_SRR7890936_50pc.vc_metrics.csv

    VARIANT CALLER SUMMARY,,Number of samples,1
    VARIANT CALLER SUMMARY,,Reads Processed,2721782043
    VARIANT CALLER SUMMARY,,Child Sample,NA
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Total,170013,100.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Biallelic,170013,100.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Multiallelic,0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNPs,138978,81.75
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Insertions (Hom),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Insertions (Het),14291,8.41
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Deletions (Hom),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Deletions (Het),16744,9.85
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Indels (Het),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Chr X number of SNPs over genome,32487
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Chr Y number of SNPs over genome,131
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,(Chr X SNPs)/(chr Y SNPs) ratio over genome,247.99
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNP Transitions,79948
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNP Transversions,59025
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Ti/Tv ratio,1.35
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Heterozygous,170013
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Homozygous,0
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Het/Hom ratio,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,In dbSNP,0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Not in dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Total,123219,100.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Biallelic,123219,100.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Multiallelic,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNPs,104900,85.13
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Insertions (Hom),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Insertions (Het),8060,6.54
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Deletions (Hom),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Deletions (Het),10259,8.33
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Indels (Het),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Chr X number of SNPs over genome,28162
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Chr Y number of SNPs over genome,8
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,(Chr X SNPs)/(chr Y SNPs) ratio over genome,3520.25
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNP Transitions,62111
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNP Transversions,42789
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Ti/Tv ratio,1.45
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Heterozygous,123219
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Homozygous,0
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Het/Hom ratio,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,In dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Not in dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Percent Callability,NA
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Percent Autosome Callability,NA
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


