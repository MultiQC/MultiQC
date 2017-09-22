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
        summaryFiles, indexSummaryFiles = self.find_log_files('interop/summary',filehandles=True), self.find_log_files('interop/index-summary',filehandles=True)
        self.runSummary,self.indexSummary = {}, {}

        for f in summaryFiles:
            log.debug( "Run Summary File Found: {}".format(f['fn'] ))      # Filename
            self.runSummary = self.parse_summary_csv(f['f'])

            #Create report Sections
            self.add_section (
                name = '{} - Run Metrics Summary'.format(f['s_name']),
                anchor = 'interop-runmetrics-summary',
                plot = self.run_metrics_summary_table(self.runSummary['summary'])
                )
            self.add_section (
                name = '{} - Run & Lane Metric Details'.format(f['s_name']),
                anchor = 'interop-runmetrics-details',
                plot = self.run_metrics_details_table(self.runSummary['details'])
                )

        for f in indexSummaryFiles:
            log.debug( "Index Summary File Found: {}".format(f['s_name'] ))      # Filename
            self.indexSummary = self.parse_index_summary_csv(f['f'])

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
            'Read 1': {
                        surface: {},
                        },
                        surface: {},
                        },
            'Read 2': {},
            ...
            }
        }
        '''
        metrics={'summary':{},
                 'details':{}
                 }

        header = []
        summary = {}
        details = {}
        section = None
        read = None
        for line in f:
            line = line.strip()
            #assume fixed file format
            #find summary header
            if line.startswith("Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30"):
                #set section to summary
                section = "summary"
                header = line.split(",")
                continue
            if section == "summary":
                data = line.split(",")
                #process summary
                summary[data[0]]={}
                for idx in range(1,len(data)):
                    summary[data[0]][header[idx]]=data[idx]
                if line.startswith("Total"):
                    section = None
                    log.debug("Finished summary")
                continue
            if line.startswith("Read") and (section is None or section == "details"):
                #set section to details
                section = "details"
                read = line
                details[read]=[]
                continue
            if line.startswith("Lane,Surface,Tiles,Density,Cluster") and section == "details":
                #get details header
                header = line.split(",")
                continue
            if section == "details":
                if line.startswith("Extracted: "):
                    section = "finish"
                    log.debug("Finished details")
                    continue
                data = line.split(",")
                #process summary
                linedata={}
                for idx in range(0,len(data)):
                    linedata[header[idx]]=data[idx]
                details[read].append(linedata)
                continue

        # import pprint
        # pp = pprint.PrettyPrinter(indent=4)
        # log.debug(pp.pprint(details))

        metrics['summary']=summary
        metrics['details']=details
        return metrics

    def parse_index_summary_csv(self,data):
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
            'no_beeswarm': True,
            'scale': False
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
            'description': 'The calculated error rate for cycles 1-35.'
        }
        headers['Error (75)'] = {
            'title': 'Error Rate 75 Cycles (%)',
            'description': 'The calculated error rate for cycles 1-75.'
        }
        headers['Error (100)'] = {
            'title': 'Error Rate 100 Cycles (%)',
            'description': 'The calculated error rate for cycles 1-100.'
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
            'col1_header': '',
            'no_beeswarm': True,
            'scale': False
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
            'no_beeswarm': True,
            'scale': False
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
            'no_beeswarm': True,
            'scale': False
        }
        return table.plot(data, headers, table_config)
