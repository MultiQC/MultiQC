#!/usr/bin/env python
""" MultiQC module to parse output from gfastats """
from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Import the gfastats submodules
from .stats import StatsReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, StatsReportMixin):
    """gfastats has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="gfastats",
            anchor="gfastats",
            target="gfastats",
            href="https://github.com/vgl-hub/gfastats",
            info=" is a single fast and exhaustive tool for summary statistics and simultaneous *fa* (fasta, fastq, gfa [.gz]) genome assembly file manipulation.",
            # doi="",  # Can't find a DOI
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = self.parse_gfastats_stats()
        if n["stats"] > 0:
            log.info("Found {} gfastats stats reports".format(n["stats"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
