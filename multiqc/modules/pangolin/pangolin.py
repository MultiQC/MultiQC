#! /usr/bin/env python

""" MultiQC module to parse output from Pangolin """

from __future__ import print_function
from collections import OrderedDict
import logging
import csv

from multiqc import config
from multiqc.plots import beeswarm, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Pangolin module"""

    def __init__(self):

        # Initialise the parent module
        super().__init__(
            name="Pangolin",
            anchor="pangolin",
            href="https://github.com/cov-lineages/pangolin",
            info="is a software package for assigning SARS-CoV-2 genome sequences to global lineages.",
        )

        # Find and parse the sample files
        self.pangolin_data = dict()
        for f in self.find_log_files("pangolin", filehandles=True):
            self.parse_pangolin_log(f)

        self.pangolin_general_stats_table()

        self.add_section(
            name="Run table",
            anchor="pangolin-run-table",
            description="Statistics gathered from the input pangolin files",
            plot=self.pangolin_table(),
        )

    def parse_pangolin_log(self, fh):
        contents = csv.DictReader(fh["f"])
        for row in contents:
            taxon_name = row["taxon"]
            row.pop("taxon")
            self.pangolin_data[taxon_name] = row

    def pangolin_general_stats_table(self):
        """Takes the parsed sample data and adds it to the general stats table"""

        headers = OrderedDict()
        headers["lineage"] = {
            "title": "Lineage",
            "description": "Lineage",
            "min": 0,
            "scale": "RdYlGn",
        }
        self.general_stats_addcols(self.pangolin_data, headers)

    def pangolin_table(self):
        """Creates the table of all data for the samples"""

        headers = OrderedDict()
        headers["lineage"] = {
            "title": "Lineage",
            "description": (
                "The most likely lineage assigned to a given sequence based on the inference engine used "
                "and the SARS-CoV-2 diversity designated."
            ),
        }
        headers["conflict"] = {
            "title": "Conflict",
            "description": (
                "In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity. "
                "If a sequence can fit into more than one category, the conflict score will be greater than 0 and reflect the number of categories the sequence could fit into. "
                "If the conflict score is 0, this means that within the current decision tree there is only one category that the sequence could be assigned to."
            ),
            "min": 0,
            "max": 1,
            "scale": "RdBl-rev",
        }

        headers["ambiguity_score"] = {
            "title": "Ambiguity score",
            "description": (
                "This score is a function of the quantity of missing data in a sequence. "
                "It represents the proportion of relevant sites in a sequence which were imputed to the reference values. "
                "A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed. "
                "This score only includes sites which are used by the decision tree to classify a sequence."
            ),
            "min": 0,
            "max": 1,
            "scale": "RdYlGn",
        }
        headers["scorpio_call"] = {
            "title": "Scorpio call",
            "description": "If a query is assigned a constellation by scorpio this call is output in this column",
        }

        headers["scorpio_support"] = {
            "title": "Scorpio support",
            "description": "The support score is the proportion of defining variants which have the alternative allele in the sequence.",
            "min": 0,
            "max": 1,
            "scale": "RdYlBl",
        }

        headers["scorpio_conflict"] = {
            "title": "Scorpio conflict",
            "description": (
                "The conflict score is the proportion of defining variants which have the reference allele in the sequence. "
                "Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores."
            ),
            "min": 0,
            "max": 1,
            "scale": "RdYlGn-rev",
        }
        headers["version"] = {
            "title": "Version",
            "description": "A version number that represents both the pango-designation number and the inference engine used to assign the lineage",
        }

        headers["pangolin_version"] = {
            "title": "Pangolin version",
            "description": "The version of pangolin software running.",
        }

        headers["pangoLEARN_version"] = {
            "title": "PangoLEARN version",
            "description": "The dated version of the pangoLEARN model installed.",
        }

        headers["pango_version"] = {
            "title": "Pango version",
            "description": "The version of pango-designation lineages that this assignment is based on.",
        }

        headers["status"] = {
            "title": "Status",
            "description": "Indicates whether the sequence passed the QC thresholds for minimum length and maximum N content.",
            "scale": "RdYlGn",
        }

        headers["note"] = {
            "title": "Note",
            "description": (
                "If any conflicts from the decision tree, this field will output the alternative assignments. "
                "If the sequence failed QC this field will describe why. "
                "If the sequence met the SNP thresholds for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb (Alternative, reference and ambiguous) alleles for that call."
            ),
        }

        return table.plot(self.pangolin_data, headers)
