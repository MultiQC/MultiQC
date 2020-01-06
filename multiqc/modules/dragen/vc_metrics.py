#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

from .utils import make_headers, Metric

# Initialise the logger
import logging
log = logging.getLogger(__name__)


class DragenVCMetrics(BaseMultiqcModule):
    def parse_vc_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/vc_metrics'):
            data = parse_vc_metrics_file(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return
        log.info('Found variant calling metrics for {} Dragen output prefixes'.format(len(data_by_sample)))

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, vc_table_headers = make_headers(all_metric_names, VC_METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, 'Variant calling')

        self.add_section(
            name='Variant calling',
            anchor='dragen-vc-metrics',
            description='Variant calling metrics. Metrics are reported for each sample in multi sample VCF '
                        'and gVCF files. Based on the run case, metrics are reported either as standard '
                        'VARIANT CALLER or JOINT CALLER. All metrics are reported for post-filter VCFs, '
                        'except for the "Filtered" metrics which represent how many variants were filtered out '
                        'from pre-filter VCF to generate the post-filter VCF.',
            plot=table.plot(data_by_sample, vc_table_headers, pconfig={'namespace': 'Variant calling'})
        )


VC_METRICS = [Metric(id, title, in_genstats=in_genstats, in_own_tabl=in_vctable, descr=descr, namespace='Variants')
              for (id, title, in_genstats, in_vctable, descr) in [
    # id_in_data                                              title (display name)   gen_stats  vc_table  description
    # Read stats:
    ('Total'                                                  , 'Variants'            , '#'  , '#'  , 'Total number of variants (SNPs + MNPs + INDELS).'),
    ('Reads Processed'                                        , 'VC reads'            , 'hid', 'hid', 'The number of reads used for variant calling, excluding any duplicate marked reads and reads falling outside of the target region'),
    ('Biallelic'                                              , 'Biallelic'           , 'hid', 'hid', 'Number of sites in a genome that contains two observed alleles, counting the reference as one, and therefore allowing for one variant allele'),
    ('Multiallelic'                                           , 'Multiallelic'        , 'hid', '%'  , 'Number of sites in the VCF that contain three or more observed alleles. The reference is counted as one, therefore allowing for two or more variant alleles'),
    ('SNPs'                                                   , 'SNPs'                , 'hid', '%'  , 'Number of SNPs in the variant set. A variant is counted as an SNP when the reference, allele 1, and allele2 are all length 1'),
    ('Indels'                                                 , 'Indels'              , 'hid', '%'  , 'Number of insetions and deletions in the variant set.'),
    ('Insertions'                                             , 'Ins'                 , 'hid', 'hid', ''),
    ('Deletions'                                              , 'Del'                 , 'hid', 'hid', ''),
    ('Insertions (Hom)'                                       , 'Hom ins'             , 'hid', 'hid', 'Number of variants that contains homozygous insertions'),
    ('Insertions (Het)'                                       , 'Het ins'             , 'hid', 'hid', 'Number of variants where both alleles are insertions, but not homozygous'),
    ('Deletions (Hom)'                                        , 'Hom del'             , 'hid', 'hid', 'Number of variants that contains homozygous deletions'),
    ('Deletions (Het)'                                        , 'Het del'             , 'hid', 'hid', 'Number of variants where both alleles are deletion, but not homozygous'),
    ('Indels (Het)'                                           , 'Het indel'           , 'hid', 'hid', 'Number of variants where genotypes are either [insertion+deletion], [insertion+snp] or [deletion+snp].'),
    ('DeNovo SNPs'                                            , 'DeNovo SNPs'         , 'hid', 'hid', 'Number of DeNovo marked SNPs, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.'),
    ('DeNovo INDELs'                                          , 'DeNovo indel'        , 'hid', 'hid', 'Number of DeNovo marked indels, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.'),
    ('DeNovo MNPs'                                            , 'DeNovo MNPs'         , 'hid', 'hid', 'Number of DeNovo marked MNPs, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.'),
    ('Chr X number of SNPs over genome'                       , 'ChrX SNP'            , 'hid', 'hid', 'Number of SNPs in chromosome X (or in the intersection of chromosome X with the target region). '
                                                                                                      'If there was no alignment to either chromosome X, this metric shows as NA'),
    ('Chr Y number of SNPs over genome'                       , 'ChrY SNP'            , 'hid', 'hid', 'Number of SNPs in chromosome Y (or in the intersection of chromosome Y with the target region). '
                                                                                                      'If there was no alignment to either chromosome Y, this metric shows as NA'),
    ('(Chr X SNPs)/(chr Y SNPs) ratio over genome'            , 'X/Y SNP ratio'       , 'hid', 'hid', 'Number of SNPs in chromosome X (or in the intersection of chromosome X with the target region) '
                                                                                                      'divided by the number of SNPs in chromosome Y (or in the intersection of chromosome Y with the '
                                                                                                      'target region). If there was no alignment to either chromosome X or chromosome Y, this metric '
                                                                                                      'shows as NA'),
    ('SNP Transitions'                                        , 'SNP Ti'              , 'hid', 'hid', 'Number of transitions - interchanges of two purines (A<->G) or two pyrimidines (C<->T)'),
    ('SNP Transversions'                                      , 'SNP Tv'              , 'hid', 'hid', 'Number of transversions - interchanges of purine and pyrimidine bases'),
    ('Ti/Tv ratio'                                            , 'Ti/Tv'               , 'hid', '#'  , 'Ti/Tv ratio: ratio of transitions to transitions.'),
    ('Heterozygous'                                           , 'Het'                 , 'hid', 'hid', 'Number of heterozygous variants'),
    ('Homozygous'                                             , 'Hom'                 , 'hid', 'hid', 'Number of homozygous variants'),
    ('Het/Hom ratio'                                          , 'Het/Hom'             , 'hid', '%'  , 'Heterozygous/ homozygous ratio'),
    ('In dbSNP'                                               , 'In dbSNP'            , 'hid', '%'  , 'Number of variants detected that are present in the dbsnp reference file. If no dbsnp file '
                                                                                                      'is provided via the --bsnp option, then both the In dbSNP and Novel metrics show as NA.'),
    ('Not in dbSNP'                                           , 'Novel'               , 'hid', 'hid', 'Number of all variants minus number of variants in dbSNP. If no dbsnp file '
                                                                                                      'is provided via the --bsnp option, then both the In dbSNP and Novel metrics show as NA.'),
    ('Percent Callability'                                    , 'Callability'         , 'hid', '#'  , 'Available only in germline mode with gVCF output. The percentage of non-N reference '
                                                                                                      'positions having a PASSing genotype call. Multi-allelic variants are not counted. '
                                                                                                      'Deletions are counted for all the deleted reference positions only for homozygous calls. '
                                                                                                      'Only autosomes and chromosomes X, Y and M are considered.'),
    ('Percent Autosome Callability'                           , 'Autosome callability', 'hid', 'hid', 'Available only in germline mode with gVCF output. The percentage of non-N reference '
                                                                                                      'positions having a PASSing genotype call. Multi-allelic variants are not counted. '
                                                                                                      'Deletions are counted for all the deleted reference positions only for homozygous calls. '
                                                                                                      'Only autosomes are considered (for all chromosomes, see the Callability metric).'),
    ('Filtered vars'                                          , 'Filtered'            , 'hid', '%'  , 'Number of raw variants minus the number of PASSed variants'),
    ('Filtered SNPs'                                          , 'Filt SNPs'           , 'hid', 'hid', 'Number of raw SNPs minus the number of PASSed SNPs'),
    ('Filtered indels'                                        , 'Filt indels'         , 'hid', 'hid', 'Number of raw indels minus the number of PASSed indels'),
]]


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

    f['s_name'] = re.search(r'(.*).vc_metrics.csv', f['fn']).group(1)

    prefilter_data = dict()
    postfilter_data = dict()

    for line in f['f'].splitlines():
        fields = line.split(',')
        analysis = fields[0]
        # sample = fields[1]
        metric = fields[2]
        value = fields[3]

        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass

        percentage = None
        if len(fields) > 4:  # percentage
            percentage = fields[4]
            try:
                percentage = float(percentage)
            except ValueError:
                pass

        if analysis == 'VARIANT CALLER SUMMARY':
            prefilter_data[metric] = value

        if analysis == 'VARIANT CALLER PREFILTER':
            prefilter_data[metric] = value

        if analysis == 'VARIANT CALLER POSTFILTER':
            postfilter_data[metric] = value
            if percentage is not None:
                postfilter_data[metric + ' pct'] = percentage

    # adding few more metrics: total insertions, deletions and indels numbers
    for data in [prefilter_data, postfilter_data]:
        data['Insertions'] = data['Insertions (Hom)'] + data['Insertions (Het)']
        data['Deletions']  = data['Deletions (Hom)']  + data['Deletions (Het)']
        data['Indels']     = data['Insertions']       + data['Deletions']
        if data['Total'] != 0:
            data['Insertions pct'] = data['Insertions'] / data['Total']
            data['Deletions pct']  = data['Deletions']  / data['Total']
            data['Indels pct']     = data['Indels']     / data['Total']

    data = postfilter_data
    # we are not really interested in all the details of pre-filtered variants, however
    # it would be nice to report how much we filtered out
    data['Filtered vars']     = prefilter_data['Total']  - data['Total']
    data['Filtered SNPs']     = prefilter_data['SNPs']   - data['SNPs']
    data['Filtered indels']   = prefilter_data['Indels'] - data['Indels']
    if prefilter_data['Total'] != 0:
        data['Filtered vars pct']   = data['Filtered vars']   / prefilter_data['Total']
    if prefilter_data['SNPs'] != 0:
        data['Filtered SNPs pct']   = data['Filtered SNPs']   / prefilter_data['SNPs']
    if prefilter_data['Indels'] != 0:
        data['Filtered indels pct'] = data['Filtered indels'] / prefilter_data['Indels']

    return data













