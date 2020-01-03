#!/usr/bin/env python
from __future__ import print_function

import itertools
import math
import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph, table

# Initialise the logger
import logging
log = logging.getLogger(__name__)

from .utils import make_headers, Metric


class DragenCoverageHist(BaseMultiqcModule):
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
                        'or a QC coverage region. The following criteria are used when calculating coverage:<br>'
                        '- Duplicate reads and clipped bases are ignored;<br>'
                        '- Only reads with MAPQ > min MAPQ and bases with BQ > min BQ are considered.',
            plot=table.plot(data_by_sample, own_tabl_headers)
        )



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

























