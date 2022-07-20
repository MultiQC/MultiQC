# !/usr/bin/env python

""" MultiQC module to parse output from DIAMOND """

from multiqc.modules.base_module import BaseMultiqcModule
import logging

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="DIAMOND",
            anchor="diamond",
            href="https://github.com/bbuchfink/diamond",
            info="a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.",
            doi="10.1038/s41592-021-01101-x",
        )

    log.info("Hello World!")

    self.find_log_files("diamond")

    # Find all files for diamond
    for f in self.find_log_files("diamond", filehandles=True):
        self.parse_logs(f)

    def parse_logs(self, f):
        file_content = logfile["f"]
        for l in file_content:
            if "Closing the output file" in l:
                print(next(file_content))
