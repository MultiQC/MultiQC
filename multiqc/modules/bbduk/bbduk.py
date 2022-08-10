#!/usr/bin/env python
""" Module to parse output from BBDuk """
from __future__ import print_function
from collections import OrderedDict
from multiqc.utils import config
from multiqc.plots import beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """BBDuk Module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="BBDuk",
            anchor="bbduk",
            href="https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/",
            info="""is a tool performing common data-quality-related trimming, 
			filtering, and masking operations with a kmer based approach""",
            ## One publication, but only for the merge tool:
            # doi="10.1371/journal.pone.0185056",
        )

        ## Define the main bbduk multiqc data object
        self.bbduk_data = dict()

        for f in self.find_log_files("bbduk", filehandles=True):
            self.parse_logs(f)

        self.bbduk_data = self.ignore_samples(self.bbduk_data)

        if len(self.bbduk_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.bbduk_data)))

        # Write data to file
        self.write_data_file(self.bbduk_data, "bbduk")

        # self.bbduk_general_stats()
        # self.bbduk_beeswarm_plot()

    def parse_logs(self, f):
        """Parses a BBDuk stdout saved in a file"""
