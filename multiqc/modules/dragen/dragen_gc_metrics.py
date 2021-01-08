import logging
import re
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


class DragenGcMetrics(BaseMultiqcModule):
    """Not to be confused with DragenFastqcGcMetrics"""

    NAMESPACE = 'Dragen GC Metrics'

    def add_gc_metrics_hist(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/gc_metrics'):
            data = parse_gc_metrics_file(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        hist_data = DragenGcMetrics.__get_normalized_gc_data(data_by_sample)
        smooth_points = 300
        self.add_section(
            name='GC Bias Histogram',
            anchor='gc-bias-hist',
            description="""
                
                """,
            plot=linegraph.plot(hist_data, {
                'id': 'gc-bias-hist',
                'title': 'Dragen: GC Bias Histogram',
                'ylab': 'Normalized Coverage',
                'xlab': '% GC',
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x} % GC</b>: {point.y} Normalized coverage',
                'smooth_points': smooth_points,
                'namespace': DragenGcMetrics.NAMESPACE,
            })
        )

        table_data = DragenGcMetrics.__get_summary_gc_data(data_by_sample)
        self.add_section(
            name='Mapping metrics',
            anchor='dragen-mapping-metrics',
            description="""
            Mapping metrics, similar to the metrics computed by the samtools-stats command.
            Shown on per read group level. To see per-sample level metrics, refer to the general
            stats table.
            """,
            plot=table.plot(table_data, pconfig={'namespace': DragenGcMetrics.NAMESPACE})
        )

        return data_by_sample.keys()

    @staticmethod
    def __get_normalized_gc_data(data_by_sample) -> dict:
        """Returns headers, data"""
        analysis = 'GC BIAS DETAILS'
        hist_data = {}
        for sample_name, sample_data in data_by_sample.items():
            # {Normalized coverage at GC 0: 0.8324, Normalized coverage at GC 1,0.9456}
            hist_data[sample_name] = {
                int(key.split()[-1]): item
                for key, item in sample_data[analysis].items()
                if key.startswith('Normalized coverage at GC')
            }
        return hist_data

    @staticmethod
    def __get_summary_gc_data(data_by_sample) -> dict:
        analysis = 'GC METRICS SUMMARY'
        summary_data = {}
        for sample_name, sample_data in data_by_sample.items():
            summary_data[sample_name] = {metric: stat for metric, stat in sample_data[analysis].items()}

        return summary_data


def parse_gc_metrics_file(f):
    """
    sample.fragment_length_hist.csv

    GC BIAS DETAILS,,Windows at GC 0,20,0.001
    GC BIAS DETAILS,,Windows at GC 1,4,0.000
    ...
    GC BIAS DETAILS,,Normalized coverage at GC 0,0.8324
    GC BIAS DETAILS,,Normalized coverage at GC 1,0.8540
    GC BIAS DETAILS,,Normalized coverage at GC 2,0.8512
    ...
    GC METRICS SUMMARY,,Window size,100
    GC METRICS SUMMARY,,Number of valid windows,3118404
    GC METRICS SUMMARY,,Number of discarded windows,366533633
    GC METRICS SUMMARY,,Average reference GC,40.84
    GC METRICS SUMMARY,,Mean global coverage,30.15
    GC METRICS SUMMARY,,Normalized coverage at GCs 0-19,1.06
    GC METRICS SUMMARY,,Normalized coverage at GCs 20-39,1.05
    GC METRICS SUMMARY,,Normalized coverage at GCs 40-59,0.97
    GC METRICS SUMMARY,,Normalized coverage at GCs 60-79,0.84
    GC METRICS SUMMARY,,Normalized coverage at GCs 80-100,0.47
    GC METRICS SUMMARY,,AT Dropout,0.58
    GC METRICS SUMMARY,,GC Dropout,2.01
    """

    f['s_name'] = re.search(r'(.*).gc_metrics.csv', f['fn']).group(1)

    data = defaultdict(dict)
    for line in f['f'].splitlines():
        tokens = line.split(',')
        # Percentage is currently unused
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
        data[analysis][metric] = stat

    return data
