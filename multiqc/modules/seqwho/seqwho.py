""" MultiQC module to parse output from SeqWho """

import json
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, bargraph
from multiqc.utils import mqc_colour

# Initialise the logger
log = logging.getLogger(__name__)


LIBRARIES = ["amplicon", "bisulf", "wxs", "chip", "wgs", "dnase", "rnaseq", "atac", "mirnaseq"]
SPECIES = ["human", "mouse"]


class MultiqcModule(BaseMultiqcModule):
    """SeqWho module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="SeqWho",
            anchor="seqwho",
            href="https://daehwankimlab.github.io/seqwho/",
            info="is a tool to determine a FASTQ(A) sequencing file identity, both source protocol and species of origin.",
            # doi=""
        )

        # Find and load any SeqWho reports
        self.seqwho_data = dict()
        self.seqwho_qualdis = dict()
        self.seqwho_qualscore = dict()
        self.seqwho_readdist = dict()

        for f in self.find_log_files("seqwho", filehandles=True):
            self.parse_json(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Filter to strip out ignored sample names
        self.seqwho_data = self.ignore_samples(self.seqwho_data)
        self.seqwho_qualdis = self.ignore_samples(self.seqwho_qualdis)
        self.seqwho_qualscore = self.ignore_samples(self.seqwho_qualscore)
        self.seqwho_readdist = self.ignore_samples(self.seqwho_readdist)

        if len(self.seqwho_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.seqwho_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.seqwho_data, "multiqc_seqwho")

        # General Stats Table
        self.seqwho_general_stats_table()

        # SeqWho MLE Species Table
        self.seqwho_species_plot()

        # SeqWho MLE Sequencing Table
        self.seqwho_library_plot()

        # Quality Distribution
        self.seqwho_qualdist_plot()

        # Quality Score
        self.seqwho_qualscore_plot()

        # Read Distribution
        self.seqwho_readdist_plot()

    # Parse our nice little JSON file
    def parse_json(self, f):
        """Parse the JSON output from SeqWho and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except json.JSONDecodeError:
            log.warning(f"Could not parse SeqWho JSON: '{f['fn']}'")
            return

        # Get sample name from JSON first
        for s_name, data in parsed_json.items():
            s_name = self.clean_s_name(s_name, f["root"])
            self.add_data_source(f, s_name)
            self.seqwho_data[s_name] = {}

            call_data = data["Call"]
            self.seqwho_data[s_name]["predicted_species"] = call_data["Species"]
            self.seqwho_data[s_name]["predicted_library"] = call_data["Library"]
            self.seqwho_data[s_name]["mle"] = call_data["Maximum Likelihood Estimate"]
            self.seqwho_data[s_name]["est_read_number"] = data["Estimated Read Number"]

            # Make MLE Matrix table for subset of predicted species and library
            mle_matrix = call_data["MLE Matrix"]
            mle_columns = SPECIES + LIBRARIES
            for i in mle_columns:
                self.seqwho_data[s_name][i] = mle_matrix[i]

            # Make a Quality Score vs count  version of the data
            self.seqwho_qualdis[s_name] = OrderedDict()
            qual_dist = data["Quality Dist"]
            for i in range(len(qual_dist)):
                self.seqwho_qualdis[s_name][i] = qual_dist[i]

            # Make a Quality score vs position  version of the data
            self.seqwho_qualscore[s_name] = OrderedDict()
            qual_scores = data["Quality Scores"]
            for i in range(len(qual_scores)):
                self.seqwho_qualscore[s_name][i] = qual_scores[i]

            # Make a Read Length vs count of the data
            self.seqwho_readdist[s_name] = OrderedDict()
            read_dist = data["Read lengths"]
            for i in range(len(read_dist)):
                self.seqwho_readdist[s_name][i] = read_dist[i]

    def seqwho_general_stats_table(self):
        """Take the parsed stats from the SeqWho report and add it to the
        General Statistics table at the top of the report"""
        species_scale = mqc_colour.mqc_colour_scale("Dark2")
        species_colors = [{v: species_scale.get_colour(i, lighten=0.5)} for i, v in enumerate(SPECIES)]
        library_scale = mqc_colour.mqc_colour_scale("Accent")
        library_colors = [{v: library_scale.get_colour(i, lighten=0.5)} for i, v in enumerate(LIBRARIES)]

        headers = {
            "predicted_species": {
                "title": "Pred Species",
                "description": "Predicted species",
                "cond_formatting_colours": species_colors,
                "cond_formatting_rules": {v: [{"s_eq": v}] for v in SPECIES},
                "scale": False,
            },
            "predicted_library": {
                "title": "Pred Library",
                "description": "Predicted library type",
                "cond_formatting_colours": library_colors,
                "cond_formatting_rules": {v: [{"s_eq": v}] for v in LIBRARIES},
                "scale": False,
            },
            "mle": {
                "title": "Max Likelihood Est",
                "description": "Overall maximum likelihood estimate (overall)",
                "min": 0,
                "format": "{:,.3f}",
                "scale": "RdYlGn",
            },
            "est_read_number": {
                "title": f"{config.read_count_prefix} Est Reads",
                "description": f"Estimated read number ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "scale": "Blues",
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.seqwho_data, headers)

    def seqwho_species_plot(self):
        """
        Take the Maximum Likelihood Species data from the SeqWho report and make
        barchart
        """

        config = {
            "id": "seqwho_species_plot",
            "title": "SeqWho: Species",
            "ylab": "Maximum likelihood estimate",
            "ymax": 1.0,
            "stacking": None,
            "cpswitch": False,
            "cpswitch_counts_label": "Maximum likelihood",
            "tt_decimals": 3,
        }

        self.add_section(
            name="Species",
            anchor="seqwho_species",
            description="This plot shows the maximum likelihood that a given species matches.",
            plot=bargraph.plot(self.seqwho_data, SPECIES, config),
        )

    def seqwho_library_plot(self):
        """Take the Maximum Likelihood Library data from the SeqWho report and make a
        barplot"""

        # Config for the plot
        config = {
            "id": "seqwho_sequencing_plot",
            "title": "SeqWho: Library Type",
            "ylab": "Maximum likelihood estimate",
            "ymax": 1.0,
            "stacking": None,
            "cpswitch": False,
            "cpswitch_counts_label": "Maximum likelihood",
            "tt_decimals": 3,
        }

        self.add_section(
            name="Library Type",
            anchor="seqwho_library",
            description="This plot shows the maximum likelihood that a given library matches.",
            plot=bargraph.plot(self.seqwho_data, LIBRARIES, config),
        )

    def seqwho_qualdist_plot(self):
        """Generate the Average Sequence Quality Distribution plot"""
        pconfig = {
            "id": "seqwho_qualdis_plot",
            "title": "SeqWho: Per Sequence Quality Scores",
            "ylab": "Reads",
            "xlab": "Sequence Quality (Phred)",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>Phred {point.x}</b>: {point.y:.0f} reads",
        }

        self.add_section(
            name="Predicted Per Sequence Quality Scores",
            anchor="seqwho_qualdis",
            description="This plot shows the predicted number of reads with average quality scores.",
            plot=linegraph.plot(self.seqwho_qualdis, pconfig),
        )

    def seqwho_qualscore_plot(self):
        """Generate the Average Quality Score per position plot"""
        pconfig = {
            "id": "seqwho_qualscore_plot",
            "title": "SeqWho: Quality Score",
            "ylab": "Sequence Quality (Phred)",
            "xlab": "Position (bp)",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>Base {point.x} </b>: {point.y}",
        }

        self.add_section(
            name="Predicted Mean Quality Scores",
            anchor="seqwho_qualscore",
            description="This plot shows the predicted mean quality value across each base position in the read.",
            plot=linegraph.plot(self.seqwho_qualscore, pconfig),
        )

    def seqwho_readdist_plot(self):
        """Generate the number of reads with certain length plot"""
        pconfig = {
            "id": "seqwho_readdist_plot",
            "title": "SeqWho: Read Distribution",
            "ylab": "Reads",
            "xlab": "Length (bp)",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>Base {point.x} bp</b>: {point.y} reads",
        }

        self.add_section(
            name="Predicted Read Distribution",
            anchor="seqwho_readdist",
            description="This plot shows predicted the number of reads with certain lengths.",
            plot=linegraph.plot(self.seqwho_readdist, pconfig),
        )
