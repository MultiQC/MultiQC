#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenPloidyEstimationMetrics(BaseMultiqcModule):
    def parse_ploidy_estimation_metrics(self):
        all_data_by_sample = dict()

        for f in self.find_log_files('dragen/ploidy_estimation_metrics'):
            data_by_sample = parse_ploidy_estimation_metrics_file(f)
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
        log.info('Found ploidy estimation metrics for {} samples'.format(len(all_data_by_sample)))

        headers = OrderedDict()
        headers['Ploidy estimation'] = {
            'title': 'Ploidy',
            'description': 'Ploidy estimation (XX, XY, X0, 00, etc.)',
            'scale': 'Set1',
        }
        self.general_stats_addcols(all_data_by_sample, headers, 'Ploidy estimation')


def parse_ploidy_estimation_metrics_file(f):
    """
    T_SRR7890936_50pc.ploidy_estimation_metrics.csv

    PLOIDY ESTIMATION,,Autosomal median coverage,55.63
    PLOIDY ESTIMATION,,X median coverage,27.44
    PLOIDY ESTIMATION,,Y median coverage,0.00
    PLOIDY ESTIMATION,,X median / Autosomal median,0.49
    PLOIDY ESTIMATION,,Y median / Autosomal median,0.00
    PLOIDY ESTIMATION,,Ploidy estimation,X0
    """

    sample = f['fn'].split('.ploidy_estimation_metrics.csv')[0]

    data_by_sample = defaultdict(dict)

    for line in f['f'].splitlines():
        _, _, metric, stat = line.split(',')
        try:
            stat = float(stat)
        except ValueError:
            pass
        data_by_sample[sample][metric] = stat

    return data_by_sample


