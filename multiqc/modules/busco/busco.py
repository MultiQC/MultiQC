#!/usr/bin/env python

""" MultiQC module to parse output from BUSCO """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ BUSCO module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='BUSCO', anchor='busco',
        href="http://busco.ezlab.org/",
        info="assesses genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs.")

        # Keys and strings, used for parsing and for plot
        self.busco_keys = {
            'complete': 'Complete BUSCOs',
            'complete_single_copy': 'Complete and single-copy BUSCOs',
            'complete_duplicated': 'Complete and duplicated BUSCOs',
            'fragmented': 'Fragmented BUSCOs',
            'missing': 'Missing BUSCOs',
            'total': 'Total BUSCO groups searched'
        }

        # Find and load any BUSCO reports
        self.busco_data = dict()
        for f in self.find_log_files('busco'):
            self.parse_busco_log(f)

        # Filter to strip out ignored sample names
        self.busco_data = self.ignore_samples(self.busco_data)

        if len(self.busco_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.busco_data)))

        # Write parsed report data to a file
        self.write_data_file(self.busco_data, 'multiqc_busco')

        # One Alignment Rate Plot per lineage
        lineages = set([ self.busco_data[s_name].get('lineage_dataset') for s_name in self.busco_data.keys() ])
        for lin in lineages:
            self.add_section (
                name = 'Lineage Assessment' if lin is None else 'Lineage: {}'.format(lin),
                anchor = 'busco-lineage-{}'.format(re.sub('\W+', '_', str(lin))),
                plot = self.busco_plot(lin)
            )

    def parse_busco_log(self, f):
        parsed_data = {}
        for l in f['f'].splitlines():
            for key, string in self.busco_keys.items():
                if string in l:
                    s = l.strip().split("\t")
                    parsed_data[key] = float(s[0])
            if 'The lineage dataset is:' in l:
                s = l.replace('# The lineage dataset is: ', '').split(' (Creation date:', 1)
                parsed_data['lineage_dataset'] = str(s[0])

        if len(parsed_data) > 0:
            self.busco_data[f['s_name']] = parsed_data

    def busco_plot (self, lin):
        """ Make the HighCharts HTML for the BUSCO plot for a particular lineage """

        data = {}
        for s_name in self.busco_data:
            if self.busco_data[s_name].get('lineage_dataset') == lin:
                data[s_name] = self.busco_data[s_name]

        plot_keys = ['complete_single_copy','complete_duplicated','fragmented','missing']
        plot_cols = ['#7CB5EC', '#434348', '#F7A35C', '#FF3C50']
        keys = OrderedDict()
        for k, col in zip(plot_keys, plot_cols):
            keys[k] = {'name': self.busco_keys[k], 'color': col}

        # Config for the plot
        config = {
            'id': 'busco_plot_{}'.format(re.sub('\W+', '_', str(lin))),
            'title': 'BUSCO Assessment Results' if lin is None else 'BUSCO Assessment Results: {}'.format(lin),
            'ylab': '# BUSCOs',
            'cpswitch_counts_label': 'Number of BUSCOs'
        }

        return bargraph.plot(data, keys, config)
