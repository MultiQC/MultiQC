#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenTimeMetrics(BaseMultiqcModule):
    def parse_time_metrics(self):
        """
        RUN TIME,,Time loading reference,00:00:00.000,0.00
        RUN TIME,,Time aligning reads,01:31:29.713,5489.71
        RUN TIME,,Time sorting and marking duplicates,01:58:18.388,7098.39
        RUN TIME,,Time DRAGStr calibration,00:00:48.280,48.28
        RUN TIME,,Time saving map/align output,01:59:24.920,7164.92
        RUN TIME,,Time partial reconfiguration,00:00:03.606,3.61
        RUN TIME,,Time variant calling,01:59:21.690,7161.69
        RUN TIME,,Time partitioning,01:31:27.108,5487.11
        RUN TIME,,Total runtime,03:32:28.426,12748.43
        """

        # TODO: figure why steps don't sum up into total runtime

        all_general_data = defaultdict(dict)
        all_bargraph_data = defaultdict(dict)
        all_bargraph_steps = list()

        for f in self.find_log_files('dragen/time_metrics'):
            bargraph_data, bargraph_steps, general_data = parse_time_metrics_file(f)
            if bargraph_data:
                for sn, data in bargraph_data.items():
                    if sn in bargraph_data:
                        log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
                    self.add_data_source(f, section='stats')
                    all_bargraph_data[sn] = data

                # TODO: if steps are differrent in different runs, make separate bargraphs
                all_bargraph_steps = bargraph_steps

        # Filter to strip out ignored sample names:
        all_bargraph_data = self.ignore_samples(all_bargraph_data)
        all_general_data = self.ignore_samples(all_general_data)

        if not all_bargraph_data:
            return
        log.info('Found time metrics for {} samples'.format(len(all_bargraph_data)))

        headers = OrderedDict()
        headers['Total runtime'] = {
            'title': 'Runtime',
            'description': 'Total runtime',
            'suffix': ' h'
        }
        self.general_stats_addcols(all_general_data, headers, 'Time metrics')

        self.add_section(
            name='Time metrics',
            anchor='dragen-time-metrics',
            description='Time metrics',
            plot=bargraph.plot(all_bargraph_data, {
                step: {'name': step} for step in all_bargraph_steps
            }, {
                'id': 'dragen_time_metrics',
                'title': 'Dragen time metrics',
                'ylab': 'hours',
                'cpswitch_counts_label': 'Hours'
            })
        )


def parse_time_metrics_file(f):
    sample = f['fn'].split('.time_metrics.csv')[0]

    data_by_sample = defaultdict(dict)
    steps = []
    general_stats = defaultdict(dict)

    for line in f['f'].splitlines():
        _, _, step, time, seconds = line.split(',')
        try:
            seconds = float(seconds)
        except ValueError:
            pass
        hours = seconds / 60 / 60
        if step == 'Total runtime':
            general_stats[sample][step] = hours
        else:
            data_by_sample[sample][step] = hours
            steps.append(step)

    return data_by_sample, steps, general_stats


