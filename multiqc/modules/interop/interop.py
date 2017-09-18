from multiqc.modules.base_module import BaseMultiqcModule
import logging
import os
import csv
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='interop', anchor='interop',
        href="http://illumina.github.io/interop/index.html",
        info="The Illumina InterOp libraries are a set of common routines used for reading and writing InterOp metric files.")

        log = logging.getLogger(__name__)
        self.find_log_files('interop/summary')

        # Find all files for mymod
        for f in self.find_log_files('interop/summary'):
            #Parse files

        #Create report Sections
        # Add section for summary stats per flow cell
        self.add_section (
            name = 'Read Metrics',
            anchor = 'interop-readmetrics',
            description = 'Metrics about each sequencing read'
            plot = self.read_metrics_table()
        )

        self.add_section (
            name = 'Lane Metrics',
            anchor = 'interop-lanemetrics',
            description = 'Metrics about each lane'
            plot = self.lane_metrics_table()
        )

    def read_metrics_table(self):
        """ Return a table with overview stats for each sequencing read for a single flow cell """
        headers = OrderedDict()
        headers['level'] = {
            'title': 'Level',
            'description': ''
        }
        headers['yield'] = {
            'title': 'Yield',
            'description': ''
        }
        headers['projected_yield'] = {
            'title': 'Projected Yield',
            'description': ''
        }
        headers['aligned'] = {
            'title': 'Aligned (%)',
            'description': ''
        }
        headers['error_rate'] = {
            'title': 'Error rate (%)',
            'description': ''
        }
        headers['intensity_c1'] = {
            'title': 'Intensity Cycle 1',
            'description': ''
        }
        headers['percent_Q30'] = {
            'title': '%>=Q30',
            'description': ''
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-readmetrics-table',
            'table_title': 'Sequencing Read Statistics',
            'col1_header': '',
            'no_beeswarm': True
        }
        return table.plot(self.read_metrics, headers, table_config)

    def lane_metrics_table(self):
        """ Return a table with overview stats for each sequencing lane for a single flow cell """
        headers = OrderedDict()
        headers['lane'] = {
            'title': 'Lane',
            'description': ''
        }
        headers['read'] = {
            'title': 'Read',
            'description': ''
        }
        headers['surface'] = {
            'title': 'Surface',
            'description': ''
        }
        headers['tiles'] = {
            'title': 'Tiles',
            'description': ''
        }
        headers['density'] = {
            'title': 'Density (K / mm^2)',
            'description': ''
        }
        headers['cluster_pf'] = {
            'title': 'Cluster PF (%)',
            'description': ''
        }
        headers['phas_prephas'] = {
            'title': 'Phas/Prephas (%)',
            'description': ''
        }
        headers['reads'] = {
            'title': 'Reads',
            'description': ''
        }
        headers['reads_pf'] = {
            'title': 'Reads PF',
            'description': ''
        }
        headers['percent_Q30'] = {
            'title': '%>=Q30',
            'description': ''
        }
        headers['yield'] = {
            'title': 'Yield',
            'description': ''
        }
        headers['cycles_error'] = {
            'title': 'Cycles Error',
            'description': ''
        }
        headers['aligned'] = {
            'title': 'Aligned (%)',
            'description': ''
        }
        headers['error'] = {
            'title': 'Error Rate (%)',
            'description': ''
        }
        headers['error_35'] = {
            'title': 'Error Rate 35 Cycles (%)',
            'description': ''
        }
        headers['error_75'] = {
            'title': 'Error Rate 75 Cycles (%)',
            'description': ''
        }
        headers['error_100'] = {
            'title': 'Error Rate 100 Cycles (%)',
            'description': ''
        }
        headers['intensity_c1'] = {
            'title': 'Intensity Cycle 1',
            'description': ''
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-lanemetrics-table',
            'table_title': 'Sequencing Lane Statistics',
            'col1_header': '',
            'no_beeswarm': True
        }
        return table.plot(self.lane_metrics, headers, table_config)


