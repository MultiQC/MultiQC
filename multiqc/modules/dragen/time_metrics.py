#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenTimeMetrics(BaseMultiqcModule):
    def parse_time_metrics(self):
        # TODO: figure why steps don't sum up into total runtime

        genstats_data_by_sample = defaultdict(dict)
        bargraph_data_by_sample = defaultdict(dict)
        all_bargraph_steps = list()

        for f in self.find_log_files('dragen/time_metrics'):
            bargraph_data, bargraph_steps, genstats_data = parse_time_metrics_file(f)
            if f['s_name'] in bargraph_data:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            bargraph_data_by_sample[f['s_name']] = bargraph_data
            genstats_data_by_sample[f['s_name']] = genstats_data

            # TODO: if steps are differrent in different runs, make separate bargraphs
            all_bargraph_steps = bargraph_steps

        # Filter to strip out ignored sample names:
        bargraph_data_by_sample = self.ignore_samples(bargraph_data_by_sample)
        genstats_data_by_sample = self.ignore_samples(genstats_data_by_sample)
        if not bargraph_data_by_sample:
            return
        log.info('Found time metrics for {} samples'.format(len(bargraph_data_by_sample)))

        headers = OrderedDict()
        headers['Total runtime'] = {
            'title': 'Run time',
            'description': 'Run time metrics',
            'suffix': ' h',
            'scale': 'Greys',
            'hidden': True,
        }
        self.general_stats_addcols(genstats_data_by_sample, headers, 'Time metrics')

        self.add_section(
            name='Run time metrics',
            anchor='dragen-time-metrics',
            description='Breakdown of the run duration of each process. '
                        'Each section on the barplot corresponds to an execution stage, with X axis values showing '
                        'time (in hours) in finished after the start of the pipeline.',
            plot=bargraph.plot(bargraph_data_by_sample, {
                step: {'name': step} for step in all_bargraph_steps
            }, {
                'id': 'dragen_time_metrics',
                'title': 'Run time metrics',
                'ylab': 'hours',
                'cpswitch_counts_label': 'Hours'
            })
        )


def parse_time_metrics_file(f):
    """
    T_SRR7890936_50pc.time_metrics.csv

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

    f['s_name'] = re.search(r'(.*).time_metrics.csv', f['fn']).group(1)

    data = defaultdict(dict)
    steps = []
    gen_stats = defaultdict(dict)

    for line in f['f'].splitlines():
        _, _, step, time, seconds = line.split(',')
        try:
            seconds = float(seconds)
        except ValueError:
            pass
        hours = seconds / 60 / 60
        if step == 'Total runtime':
            gen_stats[step] = hours
        else:
            data[step] = hours
            steps.append(step)

    return data, steps, gen_stats


