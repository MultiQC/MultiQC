#!/usr/bin/env python

""" MultiQC module to parse output from sourmash """

import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

from .compare import compare

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, compare):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="sourmash",
            anchor="sourmash",
            href="https://github.com/sourmash-bio/sourmash",
            info="Quickly search, compare, and analyze genomic and metagenomic data sets.",
            doi="10.21105/joss.00027",
        )

        n = dict()
        n["compare"] = self.parse_compare()
        if n["compare"] > 0:
            log.info("Found {} compare results".format(n["compare"]))

        if sum(n.values()) == 0:
            raise UserWarning
