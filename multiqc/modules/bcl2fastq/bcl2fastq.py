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

        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_bysample = dict()
        self.source_files = dict()

    def gather(self):
        # Gather data from all json files
        for myfile in self.find_log_files('bcl2fastq'):
            self.parse_file_as_json(myfile)

        # Collect counts by lane and sample (+source_files)
        self.split_data_by_lane_and_sample()

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)

        # Return with Warning if no files are found
        if len(self.bcl2fastq_bylane) == 0 and len(self.bcl2fastq_bysample) == 0:
            log.debug("Could not find any bcl2fastq data in {}".format(config.analysis_dir))
            raise UserWarning

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(s_name=s, source=",".join(list(set(self.source_files[s]))), module='bcl2fastq', section='bcl2fastq-bysample')

        # Add sample counts to general stats table
        self.add_general_stats()
        self.write_data_file({str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()}, 'multiqc_bcl2fastq_bylane')
        self.write_data_file(self.bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        # Add section for counts by lane
        self.add_section (
            name = 'bcl2fastq by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads)',
            helptext = "Perfect index reads are those that do not have a single mismatch. All samples of a lane are combinned. Undetermined reads are treated as a third category. To avoid conflicts the runId is prepended.",
            plot = bargraph.plot({key: {"imperfect": self.bcl2fastq_bylane[key]["total"]-self.bcl2fastq_bylane[key]["perfectIndex"], "perfect": self.bcl2fastq_bylane[key]["perfectIndex"], "undetermined": self.bcl2fastq_bylane[key]["undetermined"]} for key in self.bcl2fastq_bylane.keys()})
        )

        # Add section for counts by sample
        self.add_section (
            name = 'bcl2fastq by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample (with number of perfect index reads)',
            helptext = "Perfect index reads are those that do not have a single mismatch. All samples are aggregated across lanes combinned. Undetermined reads are ignored. Undetermined reads are treated as a separate sample. To avoid conflicts the runId is prepended.",
            plot = bargraph.plot({key: {"imperfect": self.bcl2fastq_bysample[key]["total"]-self.bcl2fastq_bysample[key]["perfectIndex"], "perfect": self.bcl2fastq_bysample[key]["perfectIndex"]} for key in self.bcl2fastq_bysample.keys()})
        )

    def parse_file_as_json(self, myfile):
        content = json.loads(myfile["f"])
        runId = content["RunId"]
        if not runId in self.bcl2fastq_data:
            self.bcl2fastq_data[runId] = dict()
        run_data = self.bcl2fastq_data[runId]
        for conversionResult in content["ConversionResults"]:
            lane = conversionResult["LaneNumber"]
            if lane in run_data:
                log.debug("Duplicate runId/lane combination found! Overwriting: {}".format(self.prepend_runid(runId, lane)))
            run_data[lane] = {"total": 0, "perfectIndex": 0, "samples": dict()}
            for demuxResult in conversionResult["DemuxResults"]:
                sample = demuxResult["SampleName"]
                run_data[lane]["samples"][sample] = {"total": 0, "perfectIndex": 0, "filename": myfile['root']+"/"+myfile["fn"]}
                run_data[lane]["total"] += demuxResult["NumberReads"]
                run_data[lane]["samples"][sample]["total"] += demuxResult["NumberReads"]
                for indexMetric in demuxResult["IndexMetrics"]:
                    run_data[lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                    run_data[lane]["samples"][sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
            run_data[lane]["samples"]["undetermined"] = {"total": conversionResult["Undetermined"]["NumberReads"], "perfectIndex": 0}

    def split_data_by_lane_and_sample(self):
        for runId in self.bcl2fastq_data.keys():
            for lane in self.bcl2fastq_data[runId].keys():
                uniqLaneName = self.prepend_runid(runId, lane)
                self.bcl2fastq_bylane[uniqLaneName] = {
                "total": self.bcl2fastq_data[runId][lane]["total"],
                "perfectIndex": self.bcl2fastq_data[runId][lane]["perfectIndex"],
                "undetermined": self.bcl2fastq_data[runId][lane]["samples"]["undetermined"]["total"]
                }
                for sample in self.bcl2fastq_data[runId][lane]["samples"].keys():
                    uniqSampleName = self.prepend_runid(runId, sample)
                    if not uniqSampleName in self.bcl2fastq_bysample:
                        self.bcl2fastq_bysample[uniqSampleName] = {"total": 0, "perfectIndex": 0}
                    self.bcl2fastq_bysample[uniqSampleName]["total"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["total"]
                    self.bcl2fastq_bysample[uniqSampleName]["perfectIndex"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["perfectIndex"]
                    if sample != "undetermined":
                        if not uniqSampleName in self.source_files:
                            self.source_files[uniqSampleName] = []
                        self.source_files[uniqSampleName].append(self.bcl2fastq_data[runId][lane]["samples"][sample]["filename"])

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

    def prepend_runid(self, runId, rest):
        return str(runId)+" - "+str(rest)
