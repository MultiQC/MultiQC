#!/usr/bin/env python
""" MultiQC module to parse output from Samtools """
from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Import the Samtools submodules
from .count import CellRangerCountMixin
from .vdj import CellRangerVdjMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, CellRangerCountMixin, CellRangerVdjMixin):
    """Cellranger has 2 main modules: count and vdj.
    This module parses data directly from the web_summary.html and summarise
    data useful for QC in the main sample table as well as generate some module specific plots
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Cellranger",
            anchor="cellranger",
            href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger",
            info="Cellranger analyze single cell expression or vdj data produced by 10X Genomics.",
            doi='10.1038/ncomms14049',
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["count"] = self.parse_count_html()
        if n["count"] > 0:
            log.info("Found {} cellranger count reports".format(n["count"]))

        n["vdj"] = self.parse_vdj_html()
        if n["vdj"] > 0:
            log.info("Found {} cellranger vdj reports".format(n["vdj"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
