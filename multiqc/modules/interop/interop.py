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
        info=" - a set of common routines used for reading and writing InterOp metric files.")

        log = logging.getLogger(__name__)

        # Parse data
        self.runSummary = {}
        self.indexSummary = {}
        for f in self.find_log_files('interop/summary',filehandles=True):
            parsed_data = self.parse_summary_csv(f['f'])
            if max(len(parsed_data['summary']), len(parsed_data['details'])) > 0:
                self.runSummary[f['s_name']] = parsed_data
        for f in self.find_log_files('interop/index-summary',filehandles=True):
            parsed_data = self.parse_index_summary_csv(f['f'])
            if max(len(parsed_data['summary']), len(parsed_data['details'])) > 0:
                self.indexSummary[f['s_name']] = parsed_data

        # No samples
        if max(len(self.runSummary), len(self.indexSummary)) == 0:
            raise UserWarning

        # Create report Sections
        if len(self.runSummary) > 0:
            self.add_section (
                name = 'Read Metrics Summary',
                anchor = 'interop-runmetrics-summary',
                description = 'Summary statistics for Total read count from each run.',
                plot = self.run_metrics_summary_table(self.runSummary)
            )
            self.add_section (
                name = 'Read Metrics per Lane',
                anchor = 'interop-runmetrics-details',
                plot = self.run_metrics_details_table(self.runSummary)
            )

        if len(self.indexSummary) > 0:
            self.add_section (
                name = 'Indexing QC Metrics summary',
                anchor = 'interop-indexmetrics-summary',
                description = 'Summary metrics about each lane',
                plot = self.index_metrics_summary_table(self.indexSummary)
            )
            self.add_section (
                name = 'Indexing QC Metrics details',
                anchor = 'interop-indexmetrics-details',
                description = ' Detail Metrics about each lane',
                plot = self.index_metrics_details_table(self.indexSummary)
            )

    def parse_summary_csv(self,f):
        metrics = {
                'summary': {},
                'details': {}
                }
        header = []
        section = None
        read = None
        for line in f:
            line = line.strip()
            # assume fixed file format
            # find summary header
            if line.startswith("Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30"):
                # set section to summary
                section = "summary"
                header = line.split(",")
            elif section == "summary":
                data = line.split(",")
                # process summary
                metrics['summary'][data[0]] = {}
                for idx in range(1,len(data)):
                    try:
                        metrics['summary'][data[0]][header[idx]] = float(data[idx])
                    except ValueError:
                        metrics['summary'][data[0]][header[idx]] = data[idx]
                if line.startswith("Total"):
                    section = None
            elif line.startswith("Read") and (section is None or section == "details"):
                # set section to details
                section = "details"
                read = line
            elif line.startswith("Lane,Surface,Tiles,Density,Cluster") and section == "details":
                # get details header
                header = line.split(",")
            elif section == "details":
                if line.startswith("Extracted: "):
                    section = "finish"
                    continue
                data = line.split(",")
                # process summary
                linedata = {}
                # Check if "surface" is total (-) else skip
                if data[1] == '-':
                    metrics['details']["{} - Lane {}".format(read,data[0])] = {}
                else:
                    continue
                for idx in range(2,len(data)):
                    try:
                        if header[idx] == 'Phas/Prephas':
                            val = data[idx].split('/')
                            linedata['Phased'] = val[0]
                            linedata['Prephased'] = val[1]
                        else:
                            linedata[header[idx]] = float(data[idx])
                    except ValueError:
                        linedata[header[idx]] = re.sub(pattern='\+/-.*', repl='', string=data[idx])
                metrics['details']["Lane {} - {}".format(data[0],read)]=linedata

        return metrics

    def parse_index_summary_csv(self,f):
        metrics={'summary':{},
                 'details':{}
                 }
        header = []
        section = None
        lane = None

        summary = {}
        details = {}

        for line in f:
            line = line.strip()
            #assume fixed file format
            if line.startswith("Lane"):
                #set lane
                lane = line
                summary[lane]={}
                continue
            if line.startswith("Total Reads,PF Reads,% Read Identified (PF),CV,Min,Max"):
                header = line.split(",")
                section = "summary"
                continue
            if line.startswith("Index Number,Sample Id,Project,Index 1 (I7),Index 2 (I5),% Read Identified (PF)"):
                header = line.split(",")
                section = "details"
                continue
            if section == "summary":
                data = line.split(",")
                for idx in range(0,len(data)):
                    summary[lane][header[idx]]=data[idx]
                continue
            if section == "details":
                data = line.split(",")
                details["{} - {}".format(data[1],lane)]={}
                for idx in range(2,len(data)):
                    details["{} - {}".format(data[1],lane)][header[idx]]=data[idx]
                continue

        metrics['summary']=summary
        metrics['details']=details

        return metrics

    def run_metrics_summary_table(self, data):

        headers = OrderedDict()
        headers['Yield'] = {
            'title': 'Yield',
            'description': 'The number of bases sequenced.',
            'scale': 'PuOr'
        }
        headers['Aligned'] = {
            'title': 'Aligned (%)',
            'description': 'The percentage of the sample that aligned to the PhiX genome',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'PiYG'
        }
        headers['Error Rate'] = {
            'title': 'Error Rate (%)',
            'description': '',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd'
        }
        headers['Intensity C1'] = {
            'title': 'Intensity Cycle 1',
            'description': ''
        }
        headers['%>=Q30'] = {
            'title': '% >= Q30',
            'description': 'Percentage of reads with quality phred score of 30 or above',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-runmetrics-summary-table',
            'table_title': 'Read metrics summary',
            'col1_header': 'Run - Read',
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]['summary']:
                tdata["{} - {}".format(s_name,key)]=data[s_name]['summary'][key]

        return table.plot(tdata, headers, table_config)

    def run_metrics_details_table(self,data):
        headers = OrderedDict()
        headers['Surface'] = {
            'title': 'Surface',
            'description': ''
        }
        headers['Tiles'] = {
            'title': 'Tiles',
            'description': 'The number of tiles per lane.',
            'hidden': True
        }
        headers['Density'] = {
            'title': 'Density',
            'description': 'The density of clusters (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.',
            'hidden': True
        }
        headers['Cluster PF'] = {
            'title': 'Cluster PF (%)',
            'description': 'The percentage of clusters passing filtering, +/- 1 standard deviation.',
            'suffix': '%',
        }
        headers['Phased'] = {
            'title': 'Phased (%)',
            'description': 'The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.',
            'format': '{:.,2f}',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd'
        }
        headers['Prephased'] = {
            'title': 'Prephased (%)',
            'description': 'The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) or jumps ahead (prephasing) the current cycle within a read.',
            'format': '{:.,2f}',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd'
        }
        headers['Reads'] = {
            'title': 'Reads',
            'description': 'The number of clusters (in millions).',
            'suffix': 'M',
        }
        headers['Reads PF'] = {
            'title': 'Reads PF',
            'description': 'The number of clusters (in millions) passing filtering.',
            'suffix': 'M',
        }
        headers['Cycles Error'] = {
            'title': 'Cycles Error',
            'description': 'The number of cycles that have been error-rated using PhiX, starting at cycle 1.'
        }
        headers['Yield'] = {
            'title': 'Yield',
            'description': 'The number of bases sequenced which passed filter.',
            'scale': 'PuOr'
        }
        headers['Aligned'] = {
            'title': 'Aligned (%)',
            'description': 'The percentage that aligned to the PhiX genome.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'PiYG'
        }
        headers['Error'] = {
            'title': 'Error Rate (%)',
            'description': 'The calculated error rate, as determined by the PhiX alignment.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd'
        }
        headers['Error (35)'] = {
            'title': 'Error Rate 35 Cycles (%)',
            'description': 'The calculated error rate for cycles 1-35.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'hidden': True
        }
        headers['Error (75)'] = {
            'title': 'Error Rate 75 Cycles (%)',
            'description': 'The calculated error rate for cycles 1-75.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'hidden': True
        }
        headers['Error (100)'] = {
            'title': 'Error Rate 100 Cycles (%)',
            'description': 'The calculated error rate for cycles 1-100.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'hidden': True
        }
        headers['Intensity C1'] = {
            'title': 'Intensity Cycle 1',
            'description': 'The corresponding intensity statistic at cycle 20 as a percentage of that value at the first cycle.',
            'help': 'equation: 100%x(Intensity at cycle 20)/(Intensity at cycle 1).'
        }
        headers['%>=Q30'] = {
            'title': '%>=Q30',
            'description': 'The percentage of bases with a quality score of 30 or higher, respectively.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-runmetrics-detail-table',
            'table_title': 'Sequencing Lane Statistics',
            'col1_header': 'Run - Lane - Read',
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]['details']:
                tdata["{} - {}".format(s_name,key)]=data[s_name]['details'][key]

        return table.plot(tdata, headers, table_config)

    def index_metrics_summary_table(self,data):
        headers = OrderedDict()
        headers['Total Reads'] = {
            'title': 'Total Reads',
            'description': 'The total number of reads for this lane (in millions).',
            'modify': lambda x: float(x) / 1000000,
            'suffix': 'M',
        }
        headers['PF Reads'] = {
            'title': 'PF Reads',
            'description': 'The total number of passing filter reads for this lane (in millions).',
            'modify': lambda x: float(x) / 1000000,
            'suffix': 'M',
        }
        headers['% Read Identified (PF)'] = {
            'title': '% Reads Identified (PF)',
            'description': 'The total fraction of passing filter reads assigned to an index.',
            'suffix': '%',
        }
        headers['CV'] = {
            'title': 'CV',
            'description': 'The coefficient of variation for the number of counts across all indexes.',
            'format': '{:.,2f}',
        }
        headers['Min'] = {
            'title': 'Min',
            'description': 'The lowest representation for any index.',
            'format': '{:.,2f}',
        }
        headers['Max'] = {
            'title': 'Max',
            'description': 'The highest representation for any index.',
            'format': '{:.,2f}',
        }
        table_config = {
            'namespace': 'interop',
            'id': 'interop-indexmetrics-summary-table',
            'table_title': 'Index Read Statistics Summary',
            'col1_header': 'Run - Lane',
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]['summary']:
                tdata["{} - {}".format(s_name,key)]=data[s_name]['summary'][key]

        return table.plot(tdata, headers, table_config)

    def index_metrics_details_table(self,data):
        headers = OrderedDict()
        headers['% Read Identified (PF)'] = {
            'title': '% Reads Identified (PF)',
            'description': 'The number of reads (only includes Passing Filter reads) mapped to this index.',
            'suffix': '%',
        }
        headers['Index 1 (I7)'] = {
            'title': 'Index 1 (I7)',
            'description': 'The sequence for the first Index Read.'
        }
        headers['Index 2 (I5)'] = {
            'title': 'Index 2 (I5)',
            'description': 'The sequence for the second Index Read.'
        }

        table_config = {
            'namespace': 'interop',
            'id': 'interop-indexmetrics-details-table',
            'table_title': 'Index Read Statistics Details',
            'col1_header': 'Run - Sample - Lane',
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]['details']:
                tdata["{} - {}".format(s_name,key)]=data[s_name]['details'][key]

        return table.plot(tdata, headers, table_config)
