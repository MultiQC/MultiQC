#! /usr/bin/env python

""" MultiQC module to parse output from Pangolin """

from __future__ import print_function
from collections import OrderedDict
import logging
import csv

from multiqc.utils import mqc_colour
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
        self.lineage_colours = dict()
        for f in self.find_log_files("pangolin", filehandles=True):
            self.parse_pangolin_log(f)

        # Filter out parsed samples based on sample name
        self.pangolin_data = self.ignore_samples(self.pangolin_data)

        # Stop if we didn't find anything
        if len(self.pangolin_data) == 0:
            raise UserWarning
        log.info("Found {} samples".format(len(self.pangolin_data)))
        self.write_data_file(self.pangolin_data, "multiqc_pangolin")

        # Assign some lineage colours
        # First, remove blank / None
        self.lineage_colours.pop("", None)
        self.lineage_colours.pop("None", None)
        cols = mqc_colour.mqc_colour_scale("Dark2", 0, len(self.lineage_colours))
        for idx, k in enumerate(self.lineage_colours):
            self.lineage_colours[k] = cols.get_colour(idx)
        # Manually add back None as grey
        self.lineage_colours["None"] = "#EFEFEF"

        self.pangolin_general_stats_table()

        self.add_section(
            name="Run table",
            anchor="pangolin-run",
            description="Statistics gathered from the input pangolin files. Hover over the column headers for descriptions and click _Help_ for more in-depth documentation.",
            helptext="""
            This table shows some of the metrics parsed by Pangolin.
            Hover over the column headers to see a description of the contents. Longer help text for certain columns is shown below:

            * **Conflict**
                * In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity.
                  If a sequence can fit into more than one category, the conflict score will be greater than `0` and reflect the number of categories the sequence could fit into.
                  If the conflict score is `0`, this means that within the current decision tree there is only one category that the sequence could be assigned to.
            * **Ambiguity score**
                * This score is a function of the quantity of missing data in a sequence.
                  It represents the proportion of relevant sites in a sequence which were imputed to the reference values.
                  A score of `1` indicates that no sites were imputed, while a score of `0` indicates that more sites were imputed than were not imputed.
                  This score only includes sites which are used by the decision tree to classify a sequence.
            * **Scorpio conflict**
                * The conflict score is the proportion of defining variants which have the reference allele in the sequence.
                  Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores.
            * **Note**
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
                # Taxon names sometimes include backslashes.
                # MultiQC assumes that these are file path separators and cleans them away
                # Bit of a nasty hack is to just replace them here first.
                taxon_name = taxon_name.replace("/", "_")
                s_name = self.clean_s_name(taxon_name, f)
                if s_name in self.pangolin_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                # Avoid generic header ID that clashes with other modules
                row["qc_status"] = row.pop("status")
                self.pangolin_data[s_name] = row
                # Just save the lineage key for now - we will sort out the colours later
                self.lineage_colours[row["lineage"]] = None
                self.lineage_colours[row["scorpio_call"]] = None
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
            "bgcols": self.lineage_colours,
        }
        self.general_stats_addcols(self.pangolin_data, headers)

    def pangolin_table(self):
        """Creates the table of all data for the samples"""

        headers = OrderedDict()
        headers["lineage"] = {
            "title": "Lineage",
            "description": """
                The most likely lineage assigned to a given sequence based on the inference engine used
                and the SARS-CoV-2 diversity designated.
            """,
            "scale": False,
            "bgcols": self.lineage_colours,
        }
        headers["conflict"] = {
            "title": "Conflict",
            "description": "Conflict between categories in decision tree",
            "min": 0,
            "max": 1,
            "scale": "RdBl-rev",
        }

        headers["ambiguity_score"] = {
            "title": "Ambiguity",
            "description": "Quantity of missing data in a sequence",
            "min": 0,
            "max": 1,
            "scale": "RdYlGn",
        }
        headers["scorpio_call"] = {
            "title": "S call",
            "description": "Scorpio: If a query is assigned a constellation by scorpio this call is output in this column",
            "scale": False,
            "bgcols": self.lineage_colours,
        }

        headers["scorpio_support"] = {
            "title": "S support",
            "description": "Scorpio: The proportion of defining variants which have the alternative allele in the sequence.",
            "min": 0,
            "max": 1,
            "scale": "RdYlBl",
        }

        headers["scorpio_conflict"] = {
            "title": "S conflict",
            "description": "Scorpio: The proportion of defining variants which have the reference allele in the sequence.",
            "min": 0,
            "max": 1,
            "scale": "RdYlGn-rev",
        }
        headers["version"] = {
            "title": "Version",
            "description": "A version number that represents both the pango-designation number and the inference engine used to assign the lineage",
            "scale": False,
            "hidden": True,
        }

        headers["pangolin_version"] = {
            "title": "Pangolin version",
            "description": "The version of pangolin software running.",
            "scale": False,
            "hidden": True,
        }

        headers["pangoLEARN_version"] = {
            "title": "PangoLEARN version",
            "description": "The dated version of the pangoLEARN model installed.",
            "scale": False,
            "hidden": True,
        }

        headers["pango_version"] = {
            "title": "Pango version",
            "description": "The version of pango-designation lineages that this assignment is based on.",
            "scale": False,
            "hidden": True,
        }

        headers["qc_status"] = {
            "title": "QC Status",
            "description": "Indicates whether the sequence passed the QC thresholds for minimum length and maximum N content.",
            "scale": False,
            "modify": lambda x: "Pass" if x == "passed_qc" else x.capitalize(),
        }

        headers["note"] = {
            "title": "Note",
            "description": "Additional information from Pangolin",
            "scale": False,
        }

        # Main table config
        table_config = {
            "namespace": "Pangolin",
            "id": "pangolin_run_table",
            "table_title": "Pangolin Run details",
        }

        return table.plot(self.pangolin_data, headers, table_config)
