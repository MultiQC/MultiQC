#!/usr/bin/env python
from __future__ import print_function

import itertools
import math
import re
from collections import OrderedDict, defaultdict

from more_itertools import flatten

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph, table

# Initialise the logger
import logging
log = logging.getLogger(__name__)

from .utils import make_headers, Metric


BASES_USED_NOTICE = 'Considering only bases usable for variant calling, i.e. excluding clipped bases, bases in duplicate reads and reads with MAPQ < min MAPQ, and bases with BQ < min BQ.'


class DragenCoverageMetrics(BaseMultiqcModule):
    def parse_coverage_metrics(self):
        data_by_phenotype_by_sample = defaultdict(dict)

        for f in self.find_log_files('dragen/wgs_coverage_metrics'):
            data_by_phenotype = parse_wgs_coverage_metrics(f)
            if f['s_name'] in data_by_phenotype:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_phenotype_by_sample[f['s_name']].update(data_by_phenotype)

        # Filter to strip out ignored sample names:
        data_by_phenotype_by_sample = self.ignore_samples(data_by_phenotype_by_sample)

        # Merge tumor and normal data:
        data_by_sample = defaultdict(dict)
        for sn in data_by_phenotype_by_sample:
            for phenotype in data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == 'normal':
                    new_sn = sn + ' normal'
                data_by_sample[new_sn] = data_by_phenotype_by_sample[sn][phenotype]

        if not data_by_sample:
            return
        log.info('Found Dragen coverage metrics for {} Dragen output prefixes'.format(len(data_by_sample)))

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, own_tabl_headers = make_headers(all_metric_names, COV_METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, 'Coverage')

        self.add_section(
            name='Coverage metrics',
            anchor='dragen-cov-metrics',
            description='Coverage metrics over a region, where the region can be the genome, a target region, '
                        'or a QC coverage region. ' + BASES_USED_NOTICE,
            plot=table.plot(data_by_sample, own_tabl_headers, pconfig={'namespace': 'Coverage'})
        )


COV_METRICS = list(flatten([[
                  Metric(id, title, in_genstats=in_genstats, in_own_tabl=in_own_tabl, unit=unit, descr=descr + ' ' + BASES_USED_NOTICE, namespace='Coverage'),
                  Metric(id.replace('{}', 'genome'), title.replace('{}', 'genome'), in_genstats=in_genstats, in_own_tabl=in_own_tabl, unit=unit, descr=descr.replace('{}', 'genome') + ' ' + BASES_USED_NOTICE, namespace='Coverage'),
                  Metric(id.replace('{}', 'genome'), title.replace('{}', 'region'), in_genstats=in_genstats, in_own_tabl=in_own_tabl, unit=unit, descr=descr.replace('{}', 'region') + ' ' + BASES_USED_NOTICE, namespace='Coverage'),
              ] for (id, title, unit, in_genstats, in_own_tabl, descr) in [
    # id_in_data                                          title (display name)    unit    gen_stats  cov_table  description
    # Read stats:
    ('Aligned reads'                                      , 'Aln reads'          ,'reads' , '#'  , '#'  , 'Total number of aligned reads.'),
    ('Aligned reads in region'                            , 'Reads on trg'       ,'reads' , 'hid', '%'  , 'Number of uniquely mapped reads to region relative to the number of uniquely mapped reads to the genome. When region is the target BED, this metric is equivalent to and replaces Capture Specificity based on target region.'),
    ('Aligned bases'                                      , 'Aln bases'          ,'bases' , '#',   '#'  , 'Total number of aligned bases.'),
    ('Aligned bases in region'                            , 'Bases on trg'       ,'bases' , '%'  , '%'  , 'Number of uniquely mapped bases to the region relative to the number of uniquely mapped bases to the genome.'),
    ('Average alignment coverage over {}'                 , 'Cov'                ,'x'     , '#'  , '#'  , 'Coverage alignment over {}: number of uniquely mapped bases to {} divided by the number of sites in {}.'),
    ('Uniformity of coverage (PCT > 0.2*mean) over {}'    , '%>0.2·mean'         ,'%'     , '#'  , '#'  , 'Percentage of sites with coverage greater than 20% of the mean coverage in region. Demonstrates the uniformity of coverage. The lower the better.'),
    ('Average chr X coverage over {}'                     , 'X cov'              ,'x'     , 'hid', 'hid', 'Average chromosome X coverage over {}. Calculated as the number of bases that aligned to the chromosome X (or to the intersection of chromosome X with the target region) divided by the total number of loci in the chromosome X (or the intersection with the target region). If there is no chromosome X in the reference genome or the region does not intersect chromosome X, this metric shows as NA.'),
    ('Average chr Y coverage over {}'                     , 'Y cov'              ,'x'     , 'hid', 'hid', 'Average chromosome Y coverage over {}. Calculated as the number of bases that aligned to the chromosome Y (or to the intersection of chromosome Y with the target region) divided by the total number of loci in the chromosome Y (or the intersection with the target region). If there is no chromosome Y in the reference genome or the region does not intersect chromosome Y, this metric shows as NA.'),
    ('Average mitochondrial coverage over {}'             , 'MT cov'             ,'x'     , 'hid', 'hid', 'Average mitochondrial chromosome coverage over {}. Calculated as the number of bases that aligned to the mitochondrial chromosome (or to the intersection of mitochondrial chromosome with the target region) divided by the total number of loci in the mitochondrial chromosome (or the intersection with the target region). If there is no mitochondrial chromosome in the reference genome or the region does not intersect mitochondrial chromosome, this metric shows as NA.'),
    ('Average autosomal coverage over {}'                 , 'Avg aut cov'        ,'x'     , 'hid', 'hid', 'Average autosomal coverage over {}. Calculated as the number of bases that aligned to the autosomal loci in {} divided by the total number of loci in the autosomal loci in {}. If there is no autosomes in the reference genome, or the region does not intersect autosomes, this metric shows as NA.'),
    ('Median autosomal coverage over {}'                  , 'Med aut cov'        ,'x'     , 'hid', 'hid', 'Median alignment coverage over the autosomal loci in {}. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.'),
    ('Mean/Median autosomal coverage ratio over {}'       , 'Mean/med aut cov ratio', ''  , 'hid', 'hid', 'Mean autosomal coverage in {} divided by the median autosomal coverage in {}. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.'),
    ('PCT of genome with coverage [100x: inf)'            , '%⩾100x'             ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 100x coverage.'),
    ('PCT of genome with coverage [ 50x: inf)'            , '%⩾50x'              ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 50x coverage.'),
    ('PCT of genome with coverage [ 20x: inf)'            , '%⩾20x'              ,'%'     , '#'  , '#'  , 'Percentage of sites in region with at least 20x coverage.'),
    ('PCT of genome with coverage [ 15x: inf)'            , '%⩾15x'              ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 15x coverage.'),
    ('PCT of genome with coverage [ 10x: inf)'            , '%⩾10x'              ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 10x coverage.'),
    ('PCT of genome with coverage [  3x: inf)'            , '%⩾3x'               ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 3x coverage.'),
    ('PCT of genome with coverage [  1x: inf)'            , '%⩾1x'               ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 1x coverage.'),
    ('PCT of genome with coverage [ 50x:100x)'            , '% 50x..100x'        ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 50x but less than 100x coverage.'),
    ('PCT of genome with coverage [ 20x: 50x)'            , '% 20x..50x'         ,'%'     , 'hid', '#'  , 'Percentage of sites in region with at least 20x but less than 50x coverage.'),
    ('PCT of genome with coverage [ 15x: 20x)'            , '% 15x..20x'         ,'%'     , 'hid', 'hid', 'Percentage of sites in region with at least 15x but less than 20x coverage.'),
    ('PCT of genome with coverage [ 10x: 15x)'            , '% 10x..15x'         ,'%'     , 'hid', 'hid', 'Percentage of sites in region with at least 10x but less than 15x coverage.'),
    ('PCT of genome with coverage [  3x: 10x)'            , '% 3x..10x'          ,'%'     , 'hid', 'hid', 'Percentage of sites in region with at least 3x but less than 10x coverage.'),
    ('PCT of genome with coverage [  1x:  3x)'            , '% 1x..3x'           ,'%'     , 'hid', 'hid', 'Percentage of sites in region with at least 1x but less than 3x coverage.'),
    ('PCT of genome with coverage [  0x:  1x)'            , '% 0x'               ,'%'     , 'hid', 'hid', 'Percentage of sites in region with no coverage.'),
]]))


def parse_wgs_coverage_metrics(f):
    """
    T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
    T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv

    The coverage metrics report outputs a _coverage_metrics.csv file, which provides metrics over a region,
    where the region can be the genome, a target region, or a QC coverage region. The first column of the output
    file contains the section name COVERAGE SUMMARY and the second column is empty for all metrics.

    The following criteria are used when calculating coverage:
    - Duplicate reads and clipped bases are ignored.
    - Only reads with MAPQ > min MAPQ and bases with BQ > min BQ are considered

    COVERAGE SUMMARY,,Aligned bases,250311219292
    COVERAGE SUMMARY,,Aligned bases in genome,250311219292,100.00
    COVERAGE SUMMARY,,Average alignment coverage over genome,82.09
    COVERAGE SUMMARY,,Uniformity of coverage (PCT > 0.2*mean) over genome,95.58
    COVERAGE SUMMARY,,PCT of genome with coverage [100x: inf),22.67
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x: inf),89.33
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: inf),95.44
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: inf),95.68
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: inf),95.92
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: inf),96.44
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x: inf),96.91
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x: inf),100.00
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x:100x),66.66
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: 50x),6.10
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: 20x),0.24
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: 15x),0.24
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: 10x),0.52
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x:  3x),0.47
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x:  1x),3.09
    COVERAGE SUMMARY,,Average chr X coverage over genome,43.24
    COVERAGE SUMMARY,,Average chr Y coverage over genome,2.88
    COVERAGE SUMMARY,,Average mitochondrial coverage over genome,25675.54
    COVERAGE SUMMARY,,Average autosomal coverage over genome,85.00
    COVERAGE SUMMARY,,Median autosomal coverage over genome,83.79
    COVERAGE SUMMARY,,Mean/Median autosomal coverage ratio over genome,1.01
    COVERAGE SUMMARY,,Aligned reads,1689488809
    COVERAGE SUMMARY,,Aligned reads in genome,1689488809,100.00
    # where "genome" can be replaced with "region"

    Add metrics into gen stats table, plus add a dedicated table like VC
    """

    data = dict()

    for line in f['f'].splitlines():
        fields = line.split(',')
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

        data[metric] = value
        if percentage is not None:
            data[metric + ' pct'] = percentage

    m = re.search(r'(.*).wgs_coverage_metrics_(\S*).csv', f['fn'])
    sample, phenotype = m.group(1), m.group(2)
    f['s_name'] = sample
    return {phenotype: data}
























