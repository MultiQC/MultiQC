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

        self.bcl2fastq_data = dict()

        for myfile in self.find_log_files('bcl2fastq'):
            content = json.loads(myfile["f"])
            for conversionResult in content["ConversionResults"]:
                for demuxResult in conversionResult["DemuxResults"]:
                    self.bcl2fastq_data[demuxResult["SampleName"]] = demuxResult["NumberReads"]

        log.info(self.bcl2fastq_data)

        if len(self.bcl2fastq_data) == 0:
            log.debug("Could not find any bcl2fastq data in {}".format(config.analysis_dir))
            raise UserWarning
