#!/usr/bin/env python
from __future__ import print_function

import itertools
import math
import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

# Initialise the logger
import logging
log = logging.getLogger(__name__)


class DragenCoverage(BaseMultiqcModule):
    def parse_coverage(self):
        perchrom_data_by_phenotype_by_sample = defaultdict(dict)

        for f in self.find_log_files('dragen/wgs_contig_mean_cov'):
            perchrom_data_by_phenotype = parse_wgs_contig_mean_cov(f)
            if f['s_name'] in perchrom_data_by_phenotype_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            perchrom_data_by_phenotype_by_sample[f['s_name']] = perchrom_data_by_phenotype

        # Filter to strip out ignored sample names:
        perchrom_data_by_phenotype_by_sample = self.ignore_samples(perchrom_data_by_phenotype_by_sample)

        # Merge tumor and normal data:
        perchrom_data_by_sample = defaultdict(dict)
        for sn in perchrom_data_by_phenotype_by_sample:
            for phenotype in perchrom_data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == 'normal':
                    new_sn = sn + ' normal'
                perchrom_data_by_sample[new_sn] = perchrom_data_by_phenotype_by_sample[sn][phenotype]

        if not perchrom_data_by_sample:
            return
        log.info('Found Dragen per-contig coverage histogram for {} Dragen output prefixes'.format(
            len(perchrom_data_by_sample)))

        self.add_section(
            name='Average coverage per contig',
            anchor='dragen-coverage-per-contig',
            description='Average coverage per contig or chromosome. Calculated as the number of bases (excluding '
                        'duplicate marked reads, reads with MAPQ=0, and clipped bases), divided by '
                        'the length of the contig or (if a target bed is used) the total length of the target '
                        'region spanning that contig>',
            plot=linegraph.plot(perchrom_data_by_sample, pconfig={
                'id': 'dragen-coverage-per-contig',
                'title': 'Average coverage per contig or chromosome',
                'ylab': 'Average coverage',
                'xlab': 'Region',
                'categories': True,
                'tt_label': '<b>{point.x}</b>: {point.y:.1f}x',
            })
        )


def parse_wgs_contig_mean_cov(f):
    """
    The Contig Mean Coverage report generates a _contig_mean_cov.csv file, which contains the estimated coverage for
    all contigs, and an autosomal estimated coverage. The file includes the following three columns

    1. Contig name
    2. Number of bases aligned to that contig, which excludes bases from duplicate marked reads, reads with MAPQ=0,
       and clipped bases.
    3. Estimated coverage, as follows: <number of bases aligned to the contig (ie, Col2)> divided by <length of the
       contig or (if a target bed is used) the total length of the target region spanning that contig>.

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

    A histogram or a plot like in mosdepth, with each chrom in X axis
    """

    perchrom_data = dict()

    for line in f['f'].splitlines():
        chrom, bases, depth = line.split(',')
        chrom = chrom.strip()
        depth = float(depth)
        # skipping unplaced and alternative contigs, as well as the mitochondria (might attract 100 times more coverage
        # than human chromosomes):
        if chrom.startswith('chrUn_') or chrom.endswith('_random') or chrom.endswith('_alt') \
                or chrom == 'chrM' or chrom == 'MT':
            continue
        perchrom_data[chrom] = depth

    def chrom_order(chrom):
        if chrom == 'Autosomal regions':
            # "Autosomal regions" average coverage goes right after all the autosomal chromosomes
            return 0
        try:
            # autosomal chromosomes go first, thus getting a negative order
            return int(chrom.replace('chr', '')) - len(perchrom_data)
        except ValueError:
            # sex and other chromosomes go in the end
            return 1

    perchrom_data = OrderedDict(sorted(perchrom_data.items(), key=lambda key_val: chrom_order(key_val[0])))

    m = re.search(r'(.*).wgs_contig_mean_cov_(\S*).csv', f['fn'])
    sample, phenotype = m.group(1), m.group(2)
    f['s_name'] = sample
    return {phenotype: perchrom_data}


def parse_wgs_coverage_metrics(f):
    """
    T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
    T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv

    COVERAGE SUMMARY,,Aligned bases,139064098059
    COVERAGE SUMMARY,,Aligned bases in genome,139064098059,100.00
    COVERAGE SUMMARY,,Average alignment coverage over genome,45.61
    COVERAGE SUMMARY,,Uniformity of coverage (PCT > 0.2*mean) over genome,95.31
    COVERAGE SUMMARY,,PCT of genome with coverage [100x: inf),0.22
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x: inf),44.01
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: inf),92.79
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: inf),94.51
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: inf),95.31
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: inf),96.09
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x: inf),96.58
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x: inf),100.00
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x:100x),43.78
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: 50x),48.79
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: 20x),1.72
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: 15x),0.80
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: 10x),0.78
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x:  3x),0.49
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x:  1x),3.42
    COVERAGE SUMMARY,,Average chr X coverage over genome,23.18
    COVERAGE SUMMARY,,Average chr Y coverage over genome,1.60
    COVERAGE SUMMARY,,Average mitochondrial coverage over genome,18694.79
    COVERAGE SUMMARY,,Average autosomal coverage over genome,47.50
    COVERAGE SUMMARY,,Median autosomal coverage over genome,48.24
    COVERAGE SUMMARY,,Mean/Median autosomal coverage ratio over genome,0.98
    COVERAGE SUMMARY,,Aligned reads,938046842
    COVERAGE SUMMARY,,Aligned reads in genome,938046842,100.00

    The coverage metrics report outputs a _coverage_metrics.csv file, which provides metrics over a region,
    where the region can be the genome, a target region, or a QC coverage region. The first column of the output
    file contains the section name COVERAGE SUMMARY and the second column is empty for all metrics.
    """

    # metrics into gen stats table, own table, plus histogram for "PCT of genome with coverage"


def parse_wgs_fine_hist(f):
    """
    T_SRR7890936_50pc.wgs_fine_hist_normal.csv
    T_SRR7890936_50pc.wgs_fine_hist_tumor.csv

    Depth,Overall
    0,104231614
    1,9430586
    2,5546235
    ...
    998,208
    999,177
    1000+,201801

    """

def parse_wgs_hist(f):
    """
    T_SRR7890936_50pc.wgs_hist_normal.csv
    T_SRR7890936_50pc.wgs_hist_tumor.csv

    PCT of bases in wgs with coverage [100x:inf), 0.22
    PCT of bases in wgs with coverage [50x:100x), 43.78
    PCT of bases in wgs with coverage [20x:50x), 48.79
    PCT of bases in wgs with coverage [15x:20x), 1.72
    PCT of bases in wgs with coverage [10x:15x), 0.80
    PCT of bases in wgs with coverage [3x:10x), 0.78
    PCT of bases in wgs with coverage [1x:3x), 0.49
    PCT of bases in wgs with coverage [0x:1x), 3.42
    """

def parse_wgs_overall_mean_cov(f):
    """
    T_SRR7890936_50pc.wgs_overall_mean_cov_normal.csv

    NORMAL Average alignment coverage over wgs, 45.61

    T_SRR7890936_50pc.wgs_overall_mean_cov_tumor.csv

    TUMOR Average alignment coverage over wgs, 82.09
    """

    # more accurate cov value. replace the value from mapping_metrics








