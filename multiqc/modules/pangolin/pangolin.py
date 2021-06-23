#! /usr/bin/env python

""" MultiQC module to parse output from Pangolin """

from __future__ import print_function
from collections import OrderedDict
import logging
import csv

from multiqc.plots import table
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
            info="uses variant calls to assign SARS-CoV-2 genome sequences to global lineages.",
        )

        # Find and parse the sample files
        self.pangolin_data = dict()
        for f in self.find_log_files("pangolin", filehandles=True):
            self.parse_pangolin_log(f)

        # Filter out parsed samples based on sample name
        self.pangolin_data = self.ignore_samples(self.pangolin_data)

        # Stop if we didn't find anything
        if len(self.pangolin_data) == 0:
            raise UserWarning
        log.info("Found {} samples".format(len(self.pangolin_data)))

        self.pangolin_general_stats_table()

        self.add_section(
            name="Run table",
            anchor="pangolin-run",
            description="Statistics gathered from the input pangolin files",
            helptext="""
            This table shows some of the metrics parsed by Pangolin.

            Hover over the column headers to see a description of the contents. Longer help text for certain columns is shown below:

            * Conflict
                * In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity.
                  If a sequence can fit into more than one category, the conflict score will be greater than 0 and reflect the number of categories the sequence could fit into.
                  If the conflict score is `0`, this means that within the current decision tree there is only one category that the sequence could be assigned to.
            * Ambiguity score
                * This score is a function of the quantity of missing data in a sequence.
                  It represents the proportion of relevant sites in a sequence which were imputed to the reference values.
                  A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed.
                  This score only includes sites which are used by the decision tree to classify a sequence.
            * Scorpio conflict
                * The conflict score is the proportion of defining variants which have the reference allele in the sequence.
                  Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores.
            * Note
                * If any conflicts from the decision tree, this field will output the alternative assignments.
                  If the sequence failed QC this field will describe why.
                  If the sequence met the SNP thresholds for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb (Alternative, reference and ambiguous) alleles for that call.
            """,
            plot=self.pangolin_table(),
        )

    def parse_pangolin_log(self, f):
        for row in csv.DictReader(f["f"]):
            try:
                taxon_name = row["taxon"]
                row.pop("taxon")
                s_name = self.clean_s_name(taxon_name, f["root"])
                if s_name in self.pangolin_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.pangolin_data[s_name] = row
            except KeyError:
                log.debug("File '{}' could not be parsed - no taxon field found.".format(f["fn"]))

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
            "description": "Conflict between categories in decision tree",
            "min": 0,
            "max": 1,
            "scale": "RdBl-rev",
        }

        headers["ambiguity_score"] = {
            "title": "Ambiguity score",
            "description": "Auantity of missing data in a sequence",
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
            "description": "The proportion of defining variants which have the alternative allele in the sequence.",
            "min": 0,
            "max": 1,
            "scale": "RdYlBl",
        }

        headers["scorpio_conflict"] = {
            "title": "Scorpio conflict",
            "description": "The proportion of defining variants which have the reference allele in the sequence.",
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
            "description": "Additional information from Pangolin",
        }

        table_config = {
            "namespace": "Pangolin",
            "id": "pangolin_run_table",
            "table_title": "Pangolin: Run details",
        }

        return table.plot(self.pangolin_data, headers, table_config)
