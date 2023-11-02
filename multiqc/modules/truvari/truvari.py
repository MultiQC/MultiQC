#!/usr/bin/env python
""" MultiQC module to parse output from truvari """
from __future__ import print_function

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# Import the truvari submodules
from .bench import BenchSummary

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, BenchSummary):
    """This is the MultiQC module for truvari."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Truvari",
            anchor="truvari",
            href="https://github.com/ACEnglish/truvari",
            info="is a toolkit for benchmarking, merging, and annotating structural variants",
            doi="https://doi.org/10.1101/2022.02.21.481353",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["bench"] = self.parse_bench_stats()
        if n["bench"] > 0:
            log.info("Found {} truvari bench reports".format(n["bench"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
