from multiqc.modules.base_module import BaseMultiqcModule
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Porechop",
            anchor="porechop",
            href="https://github.com/rrwick/Porechop",
            info="A tool for finding and removing adapters from Oxford Nanopore reads.",
        )

        log.info("Hello World!")

        self.find_log_files("porechop")

        # Find all files for porechop
        for f in self.find_log_files("porechop", filehandles=True):
            self.parse_logs(f)

    def parse_logs(self, logfile):

        file_content = logfile["f"]
        for l in file_content:
            if "Loading reads" in l:
                print(next(file_content))

        ## SAMPLE NAME
        ## Loading reads
        ## test.fastq.gz

        ## DETECTION
        ## Looking for known adapter sets
        ## 100 / 100 (100.0%)

        ## TRIMMING
        ## 100 / 100 reads had adapters trimmed from their start (6,986 bp removed)
        ## 100 / 100 reads had adapters trimmed from their end (5,389 bp removed)

        ## MIDDLE ADAPTERS
        ## Splitting reads containing middle adapters
        ## 100 / 100 (100.0%)

        ## 0 / 100 reads were split based on middle adapters
