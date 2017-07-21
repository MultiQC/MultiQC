from multiqc.modules.base_module import BaseMultiqcModule
import logging
import json
from multiqc import config

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcl2fastq', anchor='bcl2fastq',
        href="https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-18.html",
        info="bcl2fastq can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.")

        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_bysample = dict()

        for myfile in self.find_log_files('bcl2fastq'):
            content = json.loads(myfile["f"])
            for conversionResult in content["ConversionResults"]:
                lane = conversionResult["LaneNumber"]
                self.bcl2fastq_bylane[lane] = {"total": 0, "perfectIndex": 0}
                for demuxResult in conversionResult["DemuxResults"]:
                    sample = demuxResult["SampleName"]
                    if not sample in self.bcl2fastq_bysample:
                        self.bcl2fastq_bysample[sample] = {"total": 0, "perfectIndex": 0}
                    self.bcl2fastq_bylane[lane]["total"] += demuxResult["NumberReads"]
                    self.bcl2fastq_bysample[sample]["total"] += demuxResult["NumberReads"]
                    for indexMetric in demuxResult["IndexMetrics"]:
                        self.bcl2fastq_bylane[lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                        self.bcl2fastq_bysample[sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]

        log.info(self.bcl2fastq_bylane)
        log.info(self.bcl2fastq_bysample)

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)

        if len(self.bcl2fastq_bylane) == 0 and len(self.bcl2fastq_bysample) == 0:
            log.debug("Could not find any bcl2fastq data in {}".format(config.analysis_dir))
            raise UserWarning
