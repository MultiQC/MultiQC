"""MultiQC module to plot output statistics from gopeaks 

https://github.com/maxsonBraunLab/gopeaks
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import linegraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="gopeaks",
            anchor="gopeaks",
            href="https://github.com/maxsonBraunLab/gopeak",
            info="Peak caller for CUT&TAG data.",
        )

        self.macs_data = dict()
        self.data = {}
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
        self.general_stats_addcols(self.data, headers)
        
    
    def parse_gopeaks_json(self, file):
        # parse the json fields
        gopeaks_data = json.load(file['f'])
        self.data[gopeaks_data['prefix']] = gopeaks_data