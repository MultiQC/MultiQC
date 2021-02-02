import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.dragen.utils import Metric, make_headers
from multiqc.plots import table

log = logging.getLogger(__name__)

METRIC_NAMES = [
    'Unique cell-barcodes',
    'UMI threshold for passing cells',
    'Passing cells',
    'Median reads per cell',
    'Median UMIs per cell',
    'Median genes per cell',
    'Total genes detected',
]
METRICS = [
    Metric(
        id=m,
        title=' '.join([w.title() if w.islower() else w for w in m.split()]),
        in_genstats='#',
        in_own_tabl='#',
        precision=0,
    ) for m in METRIC_NAMES
]


class DragenScRnaMetrics(BaseMultiqcModule):
    def add_sc_rna_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/sc_rna_metrics'):
            data = parse_time_metrics_file(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        gen_stats_headers, table_headers = make_headers(METRIC_NAMES, METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers)
        self.add_section(
            name='Single-Cell RNA Metrics',
            anchor='sc-rna-metrics',
            description="""
            Single Cell RNA Metrics.
            """,
            plot=table.plot(
                {
                    sample_name: {
                        metric_name: int(stat) for metric_name, stat in metric.items()
                        if metric_name in METRIC_NAMES
                    }
                    for sample_name, metric in data_by_sample.items()
                },
                table_headers
            )
        )

        return data_by_sample.keys()


def parse_time_metrics_file(f):
    """
    sample.scRNA.metrics.csv

    RUN TIME,,Time loading reference,00:01:31.289,91.29
    RUN TIME,,Time aligning reads,00:00:25.190,25.19
    RUN TIME,,Time duplicate marking,00:00:01.817,1.82
    RUN TIME,,Time sorting and marking duplicates,00:00:07.368,7.37
    RUN TIME,,Time DRAGStr calibration,00:00:07.069,7.07
    """
    f['s_name'] = re.search(r'(.*).scRNA.metrics.csv', f['fn']).group(1)

    data = {}
    for line in f['f'].splitlines():
        tokens = line.split(',')
        if len(tokens) == 4:
            analysis, _, metric, stat = tokens
            percentage = None
        elif len(tokens) == 5:
            analysis, _, metric, stat, percentage = tokens
        else:
            raise ValueError(f'Unexpected number of tokens in line {line}')

        try:
            stat = float(stat)
        except ValueError:
            pass
        data[metric] = stat

    return data
