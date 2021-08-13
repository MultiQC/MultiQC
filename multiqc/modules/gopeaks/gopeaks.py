"""MultiQC module to plot output statistics from gopeaks peak caller 

https://github.com/maxsonBraunLab/gopeaks
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Gopeaks",
            anchor="Gopeaks",
            href="https://github.com/maxsonBraunLab/gopeak",
            info="is designed to call peaks from aligned CUT&TAG sequencing reads.\
                  Gopeaks uses a binomial distribution to model the read counts \
                  in sliding windows across the genome and calculate peak regions \
                  that are enriched over the background."
        )

        self.macs_data = dict()
        self.table_data = {}
        self.plot_data = {}
        for f in self.find_log_files("gopeaks", filehandles=True):
            self.add_data_source(f)
            self.parse_gopeaks_json(f)

        headers = OrderedDict({
            'gopeaks_version': {
                'title': "Gopeaks Version",
            },
            'date': {
                'title': 'Timestamp',
            },
            'elapsed': {
                'title': 'Elapsed',
            },
            'peak_counts': {
                'title': 'Peak Counts',
                'format': '{:d}',
            },
        })

        # get display config
        gopeaks_config = getattr(config, 'gopeaks_config', {})
        show_table_data = gopeaks_config.get('show_table_data', True)
        show_peak_counts_plot = gopeaks_config.get('show_peak_counts_plot', True)

        if show_table_data:
            self.general_stats_addcols(self.table_data, headers)
        if show_peak_counts_plot:
            self.add_plot_data()
        
    
    def parse_gopeaks_json(self, file):

        # parse the gopeaks metrics 
        gopeaks_data = json.load(file['f'])
        sample = gopeaks_data['prefix']

        # log the clobber
        if sample in self.table_data:
            log.debug("Duplicate sample name found, overwriting: {}.".format(sample))

        self.table_data[gopeaks_data['prefix']] = gopeaks_data
        self.plot_data[gopeaks_data['prefix']] = {"counts": gopeaks_data['peak_counts']}

    def add_plot_data(self):
        self.add_section(
            name = 'Peak Counts',
            anchor = 'gopeaks_counts',
            description = 'This plot displays the number of peaks determined by gopeaks for each sample.',
            plot = bargraph.plot(self.plot_data)
        )
