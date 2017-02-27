#!/usr/bin/env python

""" MultiQC module to parse output from BUSCO """

from __future__ import print_function
from collections import OrderedDict
import logging
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
        for f in self.find_log_files(config.sp['busco']):
            self.parse_busco_log(f)

        if len(self.busco_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.busco_data)))

        # Write parsed report data to a file
        self.write_data_file(self.busco_data, 'multiqc_busco')

        # State the lineage(s)
        lineages = set()
        for s_name in self.busco_data:
            lineages.add(self.busco_data[s_name].get('lineage_dataset'))
        if len(lineages) > 1:
            self.intro += '<p class="text-danger"><span class="glyphicon glyphicon-exclamation-sign"></span> Lineage datasets used: <em>{}</em></p>'.format(', '.join(lineages))
        else:
            self.intro += '<p class="text-info"><span class="glyphicon glyphicon-info-sign"></span> Lineage dataset used: <em>{}</em></p>'.format(lineages.pop())

        # Alignment Rate Plot
        self.intro += self.busco_plot()


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

    def busco_plot (self):
        """ Make the HighCharts HTML for the BUSCO plot """

        keys = OrderedDict()
        for k in ['complete_single_copy','complete_duplicated','fragmented','missing']:
            keys[k] = {'name': self.busco_keys[k]}

        # Config for the plot
        config = {
            'id': 'busco_plot',
            'title': 'BUSCO: Summarized benchmarking',
            'ylab': '# BUSCOs',
            'cpswitch_counts_label': 'Number of BUSCOs'
        }

        return bargraph.plot(self.busco_data, keys, config)
