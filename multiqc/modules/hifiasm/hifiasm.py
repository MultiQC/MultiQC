#!/usr/bin/env python

""" MultiQC module to parse output from HiFiasm """

import json
import logging
import re

from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="HiFiasm",
            anchor="hifiasm",
            href="https://github.com/chhylp123/hifiasm",
            info=": a haplotype-resolved assembler for accurate Hifi reads",
            doi="10.1038/s41592-020-01056-5",
        )

        # To store the mod data
        self.hifiasm_data = dict()
        self.parse_hifiasm_log_files()
        self.hifiasm_data = self.ignore_samples(self.hifiasm_data)

        # If we found no data
        if not self.hifiasm_data:
            raise UserWarning
        log.info("Found {} reports".format(len(self.hifiasm_data)))

        self.write_data_file(self.ccs_data, "multiqc_hifiasm_report")

    def parse_hifiasm_log_files(self):
        for f in self.find_log_files("hifiasm", filehandles=True):
            print(f)
