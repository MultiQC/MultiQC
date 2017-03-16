#!/usr/bin/env python

""" MultiQC module to parse output from HTSeq Count """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Adapter Removal',
        anchor='adapterRemoval', target='Adapter Removal',
        href='https://github.com/MikkelSchubert/adapterremoval',
        info=" rapid adapter trimming, identification, and read merging ")

        self.s_name = None
        self.adapter_removal_data = dict()
        self.adapter_removal_counts_mate1 = dict()
        self.adapter_removal_counts_mate2 = dict()
        self.adapter_removal_counts_singleton = dict()
        self.adapter_removal_counts_discarged = dict()
        self.adapter_removal_counts_all = dict()

        for f in self.find_log_files(config.sp['adapterRemoval'], filehandles=True):
            self.s_name = f['s_name']
            parsed_data = self.parse_settings_file(f)
            if parsed_data is not None:
                self.adapter_removal_data[self.s_name] = parsed_data

        if len(self.adapter_removal_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.adapter_removal_data)))

        # todo:
        # Write parsed report data to a file

        self.adapter_removal_stats_table()
        self.intro += self.adapter_removal_counts_chart()
        self.intro += self.adapter_removal_length_dist_plot()

    def parse_settings_file(self, f):

        self.result_data = {
            'total': None,
            'unaligned': None,
            'aligned': None,
            'percent_aligned': None,
        }

        settings_data = {'header': []}

        block_title = None
        for i, line in enumerate(f['f']):

            line = line.rstrip('\n')
            if line == '':
                continue

            if not block_title:
                block_title = 'header'
                settings_data[block_title].append(str(line))
                continue

            if line.startswith('['):
                block_title = str(line.strip('[]'))
                settings_data[block_title] = []
                continue

            settings_data[block_title].append(str(line))

        self.set_result_data(settings_data)

        # todo: return or set by member var?
        return self.result_data

    def set_result_data(self, settings_data):
        self.set_trim_stat(settings_data['Trimming statistics'])
        self.set_len_dist(settings_data['Length distribution'])

    def set_trim_stat(self, trim_data):
        data_pattern = {'total': 0, 'unaligned': 1, 'aligned': 2}

        for title, key in data_pattern.iteritems():
            tmp = trim_data[key]
            value = tmp.split(': ')[1]
            self.result_data[title] = int(value)

        self.result_data['percent_aligned'] = round((self.result_data['aligned'] * 100) / self.result_data['total'], 2)

    def set_len_dist(self, len_dist_data):
        self.adapter_removal_counts_mate1[self.s_name] = dict()
        self.adapter_removal_counts_mate2[self.s_name] = dict()
        self.adapter_removal_counts_singleton[self.s_name] = dict()
        self.adapter_removal_counts_discarged[self.s_name] = dict()
        self.adapter_removal_counts_all[self.s_name] = dict()

        for line in len_dist_data[1:]:
            #length, mate1, discarded, all = line.rstrip('\n').split('\t')
            l_data = line.rstrip('\n').split('\t')
            self.adapter_removal_counts_mate1[self.s_name][l_data[0]] = int(l_data[1])
            self.adapter_removal_counts_mate2[self.s_name][l_data[0]] = int(l_data[2])
            self.adapter_removal_counts_singleton[self.s_name][l_data[0]] = int(l_data[3])
            self.adapter_removal_counts_discarged[self.s_name][l_data[0]] = int(l_data[4])
            self.adapter_removal_counts_all[self.s_name][l_data[0]] = int(l_data[5])

    def adapter_removal_stats_table(self):

        headers = OrderedDict()
        headers['percent_aligned'] = {
            'title': '% Aligned',
            'description': '% Aligned reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
            'format': '{:.1f}%',
            'shared_key': 'percent_aligned',
        }
        headers['aligned'] = {
            'title': 'Total Aligned',
            'description': 'total aligned reads',
            'min': 0,
            'scale': 'PuBu',
            'shared_key': 'aligned'
        }
        self.general_stats_addcols(self.adapter_removal_data, headers)

    def adapter_removal_counts_chart(self):

        cats = OrderedDict()
        cats['aligned'] =      { 'name': 'aligned' }
        cats['unaligned'] =     { 'name': 'unaligned' }
        config = {
            'id': 'adapter_removal_alignment_plot',
            'title': 'Adapter Removal Alignments',
            'ylab': '# Reads',
            'hide_zero_cats': False,
            'cpswitch_counts_label': 'Number of Reads'
        }
        return bargraph.plot(self.adapter_removal_data, cats, config)

    def adapter_removal_length_dist_plot(self):

        html = '<p>add description here...</p>'

        pconfig = {
            'id': 'adapter_removal_lenght_count_plot',
            'title': 'Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Mate1', 'ylab': 'Count'},
                            {'name': 'Mate2', 'ylab': 'Count'},
                            {'name': 'Singleton', 'ylab': 'Count'},
                            {'name': 'Discarged', 'ylab': 'Count'},
                            {'name': 'All', 'ylab': 'Count'},]
        }
        #print(self.adapter_removal_counts_mate1)
        lineplot_data = [self.adapter_removal_counts_mate1,
                         self.adapter_removal_counts_mate2,
                         self.adapter_removal_counts_singleton,
                         self.adapter_removal_counts_discarged,
                         self.adapter_removal_counts_all]
        html += linegraph.plot(lineplot_data, pconfig)
        return html