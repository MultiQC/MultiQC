#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.dragen.utils import Metric
from multiqc.plots import linegraph, bargraph, table, beeswarm

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenTimeMetrics(BaseMultiqcModule):
    """
    Rendering a beeswarm. Could go with a bargraph, but the steps don't sum up to total runtime,
    meaning they are executed partially in parallel, so bars need to be overlapping.
    MultiQC doesn't support plots like that, so just showing a table.
    #TODO: figure why steps don't sum up into total runtime
    """

    def parse_time_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files('dragen/time_metrics'):
            data = parse_time_metrics_file(f)
            if f['s_name'] in data:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return
        log.info('Found time metrics for {} samples'.format(len(data_by_sample)))

        headers = OrderedDict()
        genstats_headers = OrderedDict()
        max_time = 0
        for sn, data in data_by_sample.items():
            for step, time in data.items():
                title = step
                if step.startswith('Time '):
                    title = step[5].upper() + step[6:]  # Time aligning reads -> Aligning reads
                headers[step] = {
                    'title': title,
                    'description': step,
                    'suffix': ' h',
                    'scale': 'Greys',
                }
                if step == 'Total runtime':
                    genstats_headers[step] = {
                        'title': 'Run time',
                        'description': step,
                        'suffix': ' h',
                        'scale': 'Greys',
                        'hidden': True,
                    }
                max_time = max(time, max_time)

        self.general_stats_addcols(data_by_sample, genstats_headers, 'Time metrics')

        # self.add_section(
        #     name='Run time metrics',
        #     anchor='dragen-time-metrics',
        #     description='Breakdown of the run duration of each process. '
        #                 'Each row corresponds to an execution stage, with X axis values showing '
        #                 'time (in hours) this stage took to execute. Some steps may run in parallel, '
        #                 'so the numbers don\'t nesesseraly sum up to total runtime.',
        #     plot=beeswarm.plot(data_by_sample, headers, pconfig={
        #         'namespace': 'Time metrics',
        #         'max': max_time
        #     })
        # )

        # Alternatively can consider bargraph plots. Howevew since the steps don't nesesseraly sum up to total run time
        # due to parallelization, bargraph without overlapping bars can be misleading.

        all_bargraph_steps = [
            'Total runtime',
            'Time partitioning',
            'Time loading reference'
            'Time aligning reads',
            'Time sorting and marking duplicates',
            'Time DRAGStr calibration',
            'Time saving map/align output',
            'Time partial reconfiguration',
            'Time variant calling',
        ]
        for sn, d in data_by_sample.items():
            for step in d:
                if step not in all_bargraph_steps:
                    all_bargraph_steps.append(step)
        self.add_section(
            name='Run time metrics',
            anchor='dragen-time-metrics',
            description='Breakdown of the run duration of each process. '
                        'Each section on the barplot corresponds to an execution stage, with X axis values showing '
                        'time (in hours) in finished after the start of the pipeline.',
            plot=bargraph.plot(data_by_sample, {
                step: {'name': step} for step in all_bargraph_steps
            }, {
                'id': 'dragen_time_metrics',
                'title': 'Dragen: Run time metrics',
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

    for line in f['f'].splitlines():
        _, _, step, time, seconds = line.split(',')
        try:
            seconds = float(seconds)
        except ValueError:
            pass
        hours = seconds / 60 / 60
        data[step] = hours

    return data


