from multiqc.modules.base_module import BaseMultiqcModule
import logging
import os
import csv
from collections import OrderedDict
from multiqc import config
from multiqc.plots import table
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Illumina InterOp Statistics', anchor='interop',
        href="http://illumina.github.io/interop/index.html",
        info="The Illumina InterOp libraries are a set of common routines used for reading and writing InterOp metric files.")

        log = logging.getLogger(__name__)

        # Find all files for mymod
        for f in self.find_log_files('interop/summary'):
            log.debug( "Run Summary File Found: {}".format(f['fn'] ))      # Filename
            self.parse_summary_csv(f['f'])
            self.runSummary = self.parse_summary_csv(f['f'])

        #Create report Sections
        self.add_section (
            name = 'Run Metrics Summary',
            anchor = 'interop-runmetrics-summary',
            description = 'Run metrics summary',
            plot = self.run_metrics_summary_table(self.runSummary['summary'])
        )
        self.add_section (
            name = 'Run & Lane Metrics',
            anchor = 'interop-runmetrics-details',
            description = 'Run & Lane metrics details',
            plot = self.run_metrics_details_table(self.runSummary['details'])
            )


        for f in self.find_log_files('interop/index-summary'):
            log.debug( "Index Summary File Found: {}".format(f['fn'] ))      # Filename

        self.add_section (
            name = 'Indexing QC Metrics',
            anchor = 'interop-indexmetrics-summary',
            description = 'Metrics about each lane',
            plot = self.index_metrics_summary_table(self.indexSummary['summary'])
            )
        self.add_section (
            name = 'Indexing QC Metrics',
            anchor = 'interop-indexmetrics-details',
            description = 'Metrics about each lane',
            plot = self.index_metrics_details_table(self.indexSummary['details'])
            )

    def parse_summary_csv(self,f):
        '''
        Required data structure
        data = {
            'Lane 1': {
                'var1': val1,
                'var2': val2,
            },
            'Lane 2': {
                'aligned': 1275,
                'not_aligned': 7328,
            }
        }
        '''
        metrics={'summary':{},
                 'details':{}
                 }

        return metrics
        # readMetricsDone=False
        # with open(f) as logfile:
        #     for line in logfile:
        #         #Read header, store as meta
        #         if line.startswith("# Version:"):
        #             self.meta['version'] = line.split(":")[1].strip()
        #             self.meta['run_id'] = next(logfile).strip()
        #         #Read readMetrics
        #         elif line.startswith("Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30"):
        #             header = line.strip().split(",")
        #             line=next(logfile)
        #             while line != '\n':
        #                 values=line.strip().split(",")
        #                 self.readMetrics[values[0]]={}
        #                 for idx in range(1,len(values)):
        #                     self.readMetrics[values[0]][header[idx]]=values[idx]
        #                 line = next(logfile)
        #             readMetricsDone=True
        #         #Read laneMetrics
        #         elif re.search("Read [0-3]",line.strip()) and readMetricsDone:
        #             read = line.strip()
        #             self.laneMetrics[read]={}
        #             header = next(logfile).strip().split(",")
        #             line=next(logfile)
        #             key=0
        #             while not re.search("Read [0-3]",line.strip()) and not re.search("Extracted:",line.strip()):
        #                 values=line.strip().split(",")
        #                 valuesDict={}
        #                 for idx in range(0,len(values)):
        #                     valuesDict[header[idx]]=values[idx]
        #                 self.laneMetrics[read][key]=valuesDict
        #                 key +=1
        #                 line = next(logfile)

    def parse_index-summary_csv(self,data):
        metrics={'summary':{},
                 'details':{}
                 }

        return metrics
    def run_metrics_summary_table(self,data):
        headers = OrderedDict()
        headers['Yield'] = {
            'title': 'Yield',
            'description': 'The number of bases sequenced.'
        }
        headers['Projected Yield'] = {
            'title': 'Projected Yield',
            'description': 'The projected number of bases expected to be sequenced at the end of the run.'
        }
        headers['Aligned'] = {
            'title': 'Aligned (%)',
            'description': 'The percentage of the sample that aligned to the PhiX genome, which is determined for each level or read independently.'
        }
        headers['Error Rate'] = {
            'title': 'Error rate (%)',
            'description': ''
        }
        headers['Intensity C1'] = {
            'title': 'Intensity Cycle 1',
            'description': ''
        }
        headers['%>=Q30'] = {
            'title': '%>=Q30',
            'description': ''
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-runmetrics-summary-table',
            'table_title': 'Run metrics summary',
            'col1_header': '',
            'no_beeswarm': True
        }
        return table.plot(data, headers, table_config)

    def run_metrics_details_table(self,data):
        headers = OrderedDict()
        headers['Lane'] = {
            'title': 'Lane',
            'description': ''
        }
        headers['Surface'] = {
            'title': 'Surface',
            'description': ''
        }
        headers['Tiles'] = {
            'title': 'Tiles',
            'description': 'The number of tiles per lane.'
        }
        headers['Density'] = {
            'title': 'Density (K / mm^2)',
            'description': 'The density of clusters (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.'
        }
        headers['Cluster PF'] = {
            'title': 'Cluster PF (%)',
            'description': 'The percentage of clusters passing filtering, +/- 1 standard deviation.'
        }
        headers['Phas/Prephas'] = {
            'title': 'Phas/Prephas (%)',
            'description': 'The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.'
        }
        headers['Reads'] = {
            'title': 'Reads',
            'description': 'The number of clusters (in millions).'
        }
        headers['Reads PF'] = {
            'title': 'Reads PF',
            'description': 'The number of clusters (in millions) passing filtering.'
        }
        headers['%>=Q30'] = {
            'title': '%>=Q30',
            'description': 'The percentage of bases with a quality score of 30 or higher, respectively.'
        }
        headers['Yield'] = {
            'title': 'Yield',
            'description': 'The number of bases sequenced which passed filter.'
        }
        headers['Cycles Error'] = {
            'title': 'Cycles Error',
            'description': 'The number of cycles that have been error-rated using PhiX, starting at cycle 1.'
        }
        headers['Aligned'] = {
            'title': 'Aligned (%)',
            'description': 'The percentage that aligned to the PhiX genome.'
        }
        headers['Error'] = {
            'title': 'Error Rate (%)',
            'description': 'The calculated error rate, as determined by the PhiX alignment. '
        }
        headers['Error (35)'] = {
            'title': 'Error Rate 35 Cycles (%)',
            'description': 'The calculated error rate for cycles 1–35.'
        }
        headers['Error (75)'] = {
            'title': 'Error Rate 75 Cycles (%)',
            'description': 'The calculated error rate for cycles 1–75.'
        }
        headers['Error (100)'] = {
            'title': 'Error Rate 100 Cycles (%)',
            'description': 'The calculated error rate for cycles 1–100.'
        }
        headers['Intensity C1'] = {
            'title': 'Intensity Cycle 1',
            'description': 'The corresponding intensity statistic at cycle 20 as a percentage of that value at the first cycle.',
            'help': 'equation: 100%x(Intensity at cycle 20)/(Intensity at cycle 1).'
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-runmetrics-detail-table',
            'table_title': 'Sequencing Lane Statistics',
            'col1_header': 'Run & Lane metrics details',
            'no_beeswarm': True
        }
        return table.plot(data, headers, table_config)

    def index_metrics_summary_table(self,data):
        headers = OrderedDict()
        headers['Total Reads'] = {
            'title': 'Total Reads',
            'description': 'The total number of reads for this lane.'
        }
        headers['PF Reads'] = {
            'title': 'PF Reads',
            'description': 'The total number of passing filter reads for this lane.'
        }
        headers['% Reads Identified (PF)'] = {
            'title': '% Reads Identified (PF)',
            'description': 'The total fraction of passing filter reads assigned to an index.'
        }
        headers['CV'] = {
            'title': 'The coefficient of variation for the number of counts across all indexes.',
            'description': ''
        }
        headers['Min'] = {
            'title': 'Min',
            'description': 'The lowest representation for any index.'
        }
        headers['Max'] = {
            'title': 'Max',
            'description': 'The highest representation for any index.'
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-indexmetrics-summary-table',
            'table_title': 'Index Read Statistics Summary',
            'col1_header': '',
            'no_beeswarm': True
        }
        return table.plot(data, headers, table_config)

    def index_metrics_details_table(self,data):
        headers = OrderedDict()
        headers['Sample ID'] = {
            'title': 'Sample ID',
            'description': 'The sample ID assigned to an index in the sample sheet.'
        }
        headers['Index 1 (I7)'] = {
            'title': 'Index 1 (I7)',
            'description': 'The sequence for the first Index Read.'
        }
        headers['Index 2 (I5)'] = {
            'title': 'Index 2 (I5)',
            'description': 'The sequence for the second Index Read.'
        }
        headers['% Reads Identified (PF)'] = {
            'title': '% Reads Identified (PF)',
            'description': 'The number of reads (only includes Passing Filter reads) mapped to this index.'
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-indexmetrics-details-table',
            'table_title': 'Index Read Statistics Details',
            'col1_header': '',
            'no_beeswarm': True
        }
        return table.plot(data, headers, table_config)
