from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Haplocheck",
            anchor="haplocheck",
            href="https://github.com/genepi/haplocheck/",
            info="Haplocheck detects in-sample contamination in mtDNA or WGS sequencing studies by analyzing the mitchondrial content.",
            doi=["10.1101/gr.256545.119"],
        )
        self.haplocheck_data: Dict = dict()

        for f in self.find_log_files("haplocheck", filehandles=True):
            self.parse_logs(f)

        self.haplocheck_data = self.ignore_samples(self.haplocheck_data)

        if len(self.haplocheck_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.haplocheck_data)} reports")

        # Write data to file
        self.write_data_file(self.haplocheck_data, "haplocheck")

        # trigger plots and gen stats
        ## TODO here

    def parse_logs(self, f):
        pass
