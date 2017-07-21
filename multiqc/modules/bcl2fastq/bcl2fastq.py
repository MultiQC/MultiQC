from multiqc.modules.base_module import BaseMultiqcModule
import logging
import json
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcl2fastq', anchor='bcl2fastq',
        href="https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-18.html",
        info="bcl2fastq can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.")

        self.bcl2fastq_data = dict()

        for myfile in self.find_log_files('bcl2fastq'):
            content = json.loads(myfile["f"])
            runId = content["RunId"]
            run_data = {'by_lane': dict(), 'by_sample': dict()}
            if runId in self.bcl2fastq_data:
                log.debug("Duplicate RunId found! Overwriting: {}".format(runId))
            self.bcl2fastq_data[runId] = run_data
            for conversionResult in content["ConversionResults"]:
                lane = conversionResult["LaneNumber"]
                run_data["by_lane"][lane] = {"total": 0, "perfectIndex": 0}
                for demuxResult in conversionResult["DemuxResults"]:
                    sample = demuxResult["SampleName"]
                    if not sample in run_data["by_sample"]:
                        run_data["by_sample"][sample] = {"total": 0, "perfectIndex": 0}
                    run_data["by_lane"][lane]["total"] += demuxResult["NumberReads"]
                    run_data["by_sample"][sample]["total"] += demuxResult["NumberReads"]
                    for indexMetric in demuxResult["IndexMetrics"]:
                        run_data["by_lane"][lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                        run_data["by_sample"][sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]

        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_bysample = dict()
        for runId in self.bcl2fastq_data.keys():
            for lane in self.bcl2fastq_data[runId]["by_lane"].keys():
                uniqLaneName = str(runId) + " - " + str(lane)
                if uniqLaneName in self.bcl2fastq_bylane:
                    log.debug("Duplicate lane found! Overwriting: {}".format(uniqLaneName))
                self.bcl2fastq_bylane[uniqLaneName] = self.bcl2fastq_data[runId]["by_lane"][lane]
            for sample in self.bcl2fastq_data[runId]["by_sample"].keys():
                uniqSampleName = str(runId) + " - " + str(sample)
                if uniqSampleName in self.bcl2fastq_bysample:
                    log.debug("Duplicate sample found! Overwriting: {}".format(uniqSampleName))
                self.bcl2fastq_bysample[uniqSampleName] = self.bcl2fastq_data[runId]["by_sample"][sample]

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)

        if len(self.bcl2fastq_bylane) == 0 and len(self.bcl2fastq_bysample) == 0:
            log.debug("Could not find any bcl2fastq data in {}".format(config.analysis_dir))
            raise UserWarning

        self.add_general_stats()
        self.write_data_file({str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()}, 'multiqc_bcl2fastq_bylane')
        self.write_data_file(self.bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        self.add_section (
            name = 'bcl2fastq by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads)',
#            help = "Perfect index reads are those that do not have a single mismatch. All samples of a lane are combinned. Undetermined reads are ignored.",
            plot = bargraph.plot({key: {"imperfect": self.bcl2fastq_bylane[key]["total"]-self.bcl2fastq_bylane[key]["perfectIndex"], "perfect": self.bcl2fastq_bylane[key]["perfectIndex"]} for key in self.bcl2fastq_bylane.keys()})
        )

        self.add_section (
            name = 'bcl2fastq by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample (with number of perfect index reads)',
#            help = "Perfect index reads are those that do not have a single mismatch. All samples are aggregated across lanes combinned. Undetermined reads are ignored.",
            plot = bargraph.plot({key: {"imperfect": self.bcl2fastq_bysample[key]["total"]-self.bcl2fastq_bysample[key]["perfectIndex"], "perfect": self.bcl2fastq_bysample[key]["perfectIndex"]} for key in self.bcl2fastq_bysample.keys()})
        )

    def add_general_stats(self):
        data = {key: {"total": self.bcl2fastq_bysample[key]["total"], "perfectPercent": '{0:.1f}'.format(100*self.bcl2fastq_bysample[key]["perfectIndex"]/self.bcl2fastq_bysample[key]["total"])} for key in self.bcl2fastq_bysample.keys()}
        headers = OrderedDict()
        headers['total'] = {
            'title': 'Total Reads',
            'description': 'Total number of reads for this sample as determined by bcl2fastq demultiplexing',
            'scale': 'Blues'
        }
        headers['perfectPercent'] = {
            'title': 'Perfect Index Read Percentage',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        self.general_stats_addcols(data, headers)
