#!/usr/bin/env python

""" MultiQC module to parse output from optitype """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    optitype module class, parses TSV output from optitype.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="OptiType",
            anchor="optitype",
            href="https://github.com/FRED-2/OptiType",
            info="Precision HLA typing from next-generation sequencing data ",
        )

        # Find and load any optitype reports
        self.optitype_data = dict()

        for f in self.find_log_files("optitype"):
            rows = f["f"].splitlines()
            #First col is empty / or always zero
            headers = rows[0].split("\t")[1:] 
            cols = rows[1].split("\t")[1:]

            self.optitype_data[f["s_name"]] = dict()
            for header, col in zip(headers, cols):
                self.optitype_data[f["s_name"]].update({header: col})

        if len(self.optitype_data) == 0:
            raise UserWarning

        self.optitype_data = self.ignore_samples(self.optitype_data)

        if len(self.optitype_data) == 0:
            raise UserWarning
            if parsed_data is not None and len(parsed_data) > 0:
                self.optitype_data[f["s_name"]] = parsed_data
                self.add_data_source(f, f["s_name"])

        # Filter to strip out ignored sample names
        self.optitype_data = self.ignore_samples(self.optitype_data)

        if len(self.optitype_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.optitype_data)))

        # Write parsed report data to a file
        self.write_data_file(self.optitype_data, "multiqc_optitype")

        #Create OptiType overview table
        self.addSummaryMetrics()


    def addSummaryMetrics(self):
        """ Take the parsed entries from OptiType and add them to the main plot """
        # A1, A2, B1, B2, C1, C2, Reads, Objective
        headers = OrderedDict()

        headers["A1"] = {
            "title": "HLA-A1",
            "description": "HLA-I A1 subtype",
            "scale": "RdGn",
            "shared_key": "hla",
        }
        headers["A2"] = {
            "title": "HLA-A2",
            "description": "HLA-I A2 subtype",
            "scale": "RdYlGn",
            "shared_key": "hla",
            "hidden": True,
        }
        headers["B1"] = {
            "title": "HLA-B1",
            "description": "HLA-I B1 subtype",
            "scale": "YlGn",
            "shared_key": "hla",
        }
        headers["B2"] = {
            "title": "HLA-B2",
            "description": "HLA-I B2 subtype",
            "scale": "OrRd",
            "shared_key": "hla",
            "hidden": True,
        }
        headers["C1"] = {
            "title": "HLA-C1",
            "description": "HLA-I C1 subtype",
            "scale": "OrRd",
            "shared_key": "hla",
        }
        headers["C2"] = {
            "title": "HLA-C2",
            "description": "HLA-I C2 subtype",
            "scale": "OrRd",
            "shared_key": "hla",
            "hidden": True,
        }
        headers["Reads"] = {
            "title": "Reads",
            "description": "Number of reads covering the HLA",
            "scale": "OrRd",
            "shared_key": "hla",
            "hidden": True,
            "format": "{:,.0f}",
        }
        headers["Objective"] = {
            "title": "Objective Score",
            "description": "Score of the objective function for the prediction.",
            "scale": "OrRd",
            "shared_key": "hla",
            "hidden": True,
        }
        self.general_stats_addcols(self.optitype_data, headers)