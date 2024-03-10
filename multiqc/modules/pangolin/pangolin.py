#! /usr/bin/env python

""" MultiQC module to parse output from Pangolin """


import csv
import logging
from typing import Optional

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.utils import mqc_colour

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
            doi="10.1093/ve/veab064",
        )

        # Find and parse the sample files
        self.pangolin_data = dict()
        self.lineage_colours = dict()
        for f in self.find_log_files("pangolin", filehandles=True):
            self.parse_pangolin_log(f)
            self.add_data_source(f)

        # Filter out parsed samples based on sample name
        self.pangolin_data = self.ignore_samples(self.pangolin_data)

        # Stop if we didn't find anything
        if len(self.pangolin_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.pangolin_data)} samples")
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
        def add_version_not_none(version, sample, name):
            if version is not None:
                self.add_software_version(version, sample, name)

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
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                # Avoid generic header ID that clashes with other modules
                if "qc_status" not in row:
                    row["qc_status"] = row.pop("status")
                if "qc_notes" in row:
                    row["qc_notes"] = _format_qc_notes(row["qc_notes"])

                self.pangolin_data[s_name] = row
                # Just save the lineage key for now - we will sort out the colours later
                self.lineage_colours[row["lineage"]] = None
                self.lineage_colours[row["scorpio_call"]] = None

                # Version info
                # Note: Excluded "version" field from software versions as this refers to
                # how the reference data was prepared. This info is still available
                # in the "Run table" table
                add_version_not_none(row.get("pangolin_version"), s_name, self.name)
                add_version_not_none(row.get("pango_version"), s_name, "Pango")
                add_version_not_none(row.get("pangoLEARN_version"), s_name, "PangoLEARN")
                add_version_not_none(row.get("scorpio_version"), s_name, "Scorpio")
                # constellation_version is sometimes "TRUE" or "FALSE" - ignore these
                constellation_version = row.get("constellation_version")
                if constellation_version not in {None, "TRUE", "FALSE"}:
                    self.add_software_version(constellation_version, s_name, "Constellations")
            except KeyError:
                log.debug(f"File '{f['fn']}' could not be parsed - no taxon field found.")

    def pangolin_general_stats_table(self):
        """Takes the parsed sample data and adds it to the general stats table"""

        headers = {
            "lineage": {
                "title": "Lineage",
                "description": "Lineage",
                "min": 0,
                "scale": "RdYlGn",
                "bgcols": self.lineage_colours,
            }
        }
        self.general_stats_addcols(self.pangolin_data, headers)

    def pangolin_table(self):
        """Creates the table of all data for the samples"""

        headers = {
            "lineage": {
                "title": "Lineage",
                "description": """
                The most likely lineage assigned to a given sequence based on the inference engine used
                and the SARS-CoV-2 diversity designated
                """,
                "scale": False,
                "bgcols": self.lineage_colours,
            },
            "conflict": {
                "title": "Conflict",
                "description": "Conflict between categories in decision tree",
                "min": 0,
                "max": 1,
                "scale": "RdBu-rev",
            },
            "ambiguity_score": {
                "title": "Ambiguity",
                "description": "Quantity of missing data in a sequence",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
            "scorpio_call": {
                "title": "S call",
                "description": "Scorpio: If a query is assigned a constellation by scorpio this call is output in this column",
                "scale": False,
                "bgcols": self.lineage_colours,
            },
            "scorpio_support": {
                "title": "S support",
                "description": "Scorpio: The proportion of defining variants which have the alternative allele in the sequence",
                "min": 0,
                "max": 1,
                "scale": "RdYlBu",
            },
            "scorpio_conflict": {
                "title": "S conflict",
                "description": "Scorpio: The proportion of defining variants which have the reference allele in the sequence",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn-rev",
            },
            "version": {
                "title": "Version",
                "description": "A version number that represents both the pango-designation number and the inference engine used to assign the lineage",
                "scale": False,
                "hidden": True,
            },
            "pangolin_version": {
                "title": "Pangolin version",
                "description": "The version of pangolin software running",
                "scale": False,
                "hidden": True,
            },
            "scorpio_version": {
                "title": "Scorpio version",
                "description": "The version of the scorpio software installed",
                "scale": False,
                "hidden": True,
            },
            "constellation_version": {
                "title": "Constellations version",
                "description": "The version of Constellations that scorpio has used to curate the lineage assignment",
                "scale": False,
                "hidden": True,
            },
            "qc_status": {
                "title": "QC Status",
                "description": "Indicates whether the sequence passed the QC thresholds for minimum length and maximum N content",
                "scale": False,
                "modify": lambda x: "Pass" if x == "passed_qc" else x.capitalize(),
            },
            "qc_notes": {
                "title": "QC Note",
                "description": "Notes specific to the QC checks run on the sequences",
                "scale": False,
            },
            "note": {
                "title": "Note",
                "description": "Additional information from Pangolin",
                "scale": False,
            },
        }

        # Main table config
        table_config = {
            "namespace": "Pangolin",
            "id": "pangolin_run_table",
            "table_title": "Pangolin Run details",
        }

        return table.plot(self.pangolin_data, headers, table_config)


def _format_qc_notes(raw: str) -> Optional[str]:
    """
    Parses QC notes, they appear to come from:
    https://github.com/cov-lineages/pangolin/blob/361f49cbffbf26eb28bed2f4a4c0e7f3d5a054cc/pangolin/utils/preprocessing.py#L91-L97
    https://github.com/cov-lineages/pangolin/blob/361f49cbffbf26eb28bed2f4a4c0e7f3d5a054cc/pangolin/utils/preprocessing.py#L179
    """
    if raw.startswith("Ambiguous_content"):
        # e.g. Ambiguous_content:0.03
        split = raw.split(":")
        if len(split) != 2:
            logging.warning(f"Expected label of format 'Ambiguous_content:0.01', found: '{raw}'")
            return None
        proportion_n = float(split[1])
        percent_n = int(proportion_n * 100)
        return f"Ambiguous content: {percent_n}%"

    # Unrecognized notes, just return them, capitalized
    return raw.capitalize()
