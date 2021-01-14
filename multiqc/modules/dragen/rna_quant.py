import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class DragenRnaMetrics(BaseMultiqcModule):
    def add_rna_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/rna'):
            data = parse_time_metrics_file(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        coding_fragments_metric_names = {
            'Transcript fragments',
            'Unknown transcript fragments',
            'Intron fragments',
            'Intergenic fragments'
        }

        strand_metrics = DragenRnaMetrics.__get_strand_metrics(data_by_sample)
        self.add_section(
            name='RNA Metrics',
            anchor='dragen-rna-metrics',
            description='RNA metrics for DRAGEN.',
            plot=bargraph.plot(
                [
                    {
                        sample: {
                            metric: stat for metric, stat in data.items() if
                            metric in coding_fragments_metric_names
                        } for sample, data in data_by_sample.items()
                    },
                    strand_metrics,

                ],
                pconfig={
                    'id': 'time_metrics_plot',
                    'title': 'Dragen: RNA Metrics',
                    'ylab': 'Fragments',
                    'cpswitch_counts_label': 'Fragments',
                    'data_labels': [
                        {
                            'name': 'Transcript / Intronic / Intergenic Fragments',
                            'ylab': 'Fragments',
                            'cpswitch_counts_label': 'Fragments',
                        },
                        {
                            'name': 'Fragment Orientations',
                            'ylab': 'Fragments',
                            'cpswitch_counts_label': 'Fragments',
                        },

                    ]
                })
        )

        return data_by_sample.keys()

    @staticmethod
    def __get_strand_metrics(data_by_sample: dict) -> dict:
        strand_metrics = {}
        for sample, data in data_by_sample.items():
            if 'Forward transcript fragments' in data:
                fragment_orientation_metric_names = {
                    'Forward transcript fragments',
                    'Reverse transcript fragments',
                    'Orientation filtered fragments',
                    'Ambiguous orientation fragments',
                }
            elif 'Strand mismatched fragments' in data:
                fragment_orientation_metric_names = {
                    'Transcript fragments',
                    'Strand mismatched fragments',
                }
            else:
                continue

            strand_metrics[sample] = {metric: stat for metric, stat in data.items()
                                      if metric in fragment_orientation_metric_names}

        return strand_metrics


def parse_time_metrics_file(f):
    """
    sample.time_metrics.csv

    RUN TIME,,Time loading reference,00:01:31.289,91.29
    RUN TIME,,Time aligning reads,00:00:25.190,25.19
    RUN TIME,,Time duplicate marking,00:00:01.817,1.82
    RUN TIME,,Time sorting and marking duplicates,00:00:07.368,7.37
    RUN TIME,,Time DRAGStr calibration,00:00:07.069,7.07
    """
    f['s_name'] = re.search(r'(.*).quant.metrics.csv', f['fn']).group(1)

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
