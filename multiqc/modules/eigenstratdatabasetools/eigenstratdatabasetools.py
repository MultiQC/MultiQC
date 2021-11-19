#!/usr/bin/env python

""" MultiQC module to parse output from eigenstrat_snp_coverage """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc.plots import bargraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """eigenstratdatabasetools module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="eigenstratdatabasetools",
            anchor="eigenstrat",
            href="https://github.com/TCLamnidis/EigenStratDatabaseTools",
            info="A set of tools to compare and manipulate the contents of EingenStrat databases, and to calculate SNP coverage statistics in such databases.",
            # No publication / DOI // doi=
        )

        # Find and load any DeDup reports
        self.snp_cov_data = dict()

        # Find and load JSON file
        for f in self.find_log_files("eigenstratdatabasetools", filehandles=True):
            self.parse_data(f)

        # Filter samples
        self.snp_cov_data = self.ignore_samples(self.snp_cov_data)

        # Return if no samples found
        if len(self.snp_cov_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.snp_cov_data)))

        # Save data output file
        self.write_data_file(self.snp_cov_data, "multiqc_snp_cov_metrics")

        # Add to General Statistics
        self.addSummaryMetrics()

    def parse_data(self, f):
        try:
            data = json.load(f["f"])
        except Exception as e:
            log.debug(e)
            log.warning("Could not parse eigenstrat_snp_coverage JSON: '{}'".format(f["fn"]))
            return

        # Parse JSON data to a dict
        for s_name in data:
            if s_name == "Metadata":
                continue

            s_clean = self.clean_s_name(s_name, f)
            if s_clean in self.snp_cov_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_clean))

            self.add_data_source(f, s_clean)
            self.snp_cov_data[s_clean] = dict()

            for k, v in data[s_name].items():
                try:
                    ## number of covered and total SNPs are integer values
                    self.snp_cov_data[s_clean][k] = int(v)
                except ValueError:
                    self.snp_cov_data[s_clean][k] = v

    def addSummaryMetrics(self):
        """Take the parsed stats from eigenstrat_snp_coverage and add it to the main plot"""

        headers = OrderedDict()
        headers["Covered_Snps"] = {
            "title": "Covered SNPs",
            "description": "The number of SNPs for which a genotype has been called.",
            "scale": "PuBuGn",
            "format": "{:,.0f}",
            "shared_key": "snp_call",
        }
        headers["Total_Snps"] = {
            "title": "Total SNPs",
            "description": "The total number of SNPs in the genotype dataset.",
            "scale": "PuBuGn",
            "format": "{:,.0f}",
            "hidden": True,
            "shared_key": "snp_call",
        }

        self.general_stats_addcols(self.snp_cov_data, headers)
