# !/usr/bin/env python

""" MultiQC module to parse output from SeqWho """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import json
import numpy as np

from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import report

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """SeqWho module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="SeqWho",
            anchor="seqwho",
            href="https://daehwankimlab.github.io/seqwho/",
            info="is a tool to determine a FASTQ(A) sequencing file identity, both source protocol and species of origin.",
        )

        # Find and load any SeqWho reports
        self.seqwho_data = dict()
        self.seqwho_qualdis = dict()
        self.seqwho_qualscore = dict()
        self.seqwho_readdist = dict()

        for f in self.find_log_files("seqwho", filehandles=True):
            try:
                self.parseJSON(f)
            except KeyError:
                logging.warning("Error loading file {}".format(f["fn"]))

        # Filter to strip out ignored sample names
        self.seqwho_data = self.ignore_samples(self.seqwho_data)
        self.seqwho_qualdis = self.ignore_samples(self.seqwho_qualdis)
        self.seqwho_qualscore = self.ignore_samples(self.seqwho_qualscore)
        self.seqwho_readdist = self.ignore_samples(self.seqwho_readdist)

        if len(self.seqwho_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.seqwho_data)))

        # Write parsed report data to a file
        self.write_data_file(self.seqwho_data, "multiqc_seqwho")

        # General Stats Table
        self.seqwho_general_stats_table()

        # SeqWho MLE Species Table
        self.add_section(
            description="This plot shows the Maximum Likelihood that a given Species matches.",
            plot=self.seqwho_species_plot(),
        )

        # SeqWho MLE Sequencing Table
        self.add_section(
            description="This plot shows the Maximum Likelihood that a given Library matches.",
            plot=self.seqwho_library_plot(),
        )

        # Quality Distribution
        self.seqwho_qualdist_plot()

        # Quality Score
        self.seqwho_qualscore_plot()

        # Read Distribution
        self.seqwho_readdist_plot()

    # Parse our nice little JSON file
    def parseJSON(self, f):
        """Parse the JSON output from SeqWho and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except:
            log.warning("Could not parse SeqWho JSON: '{}'".format(f["fn"]))
            return None

        # Get sample name from JSON first
        sample_names = parsed_json.keys()
        for s_name in sample_names:
            self.add_data_source(f, s_name)
            self.seqwho_data[s_name] = {}

            call_data = parsed_json[s_name]["Call"]
            self.seqwho_data[s_name]["predicted_species"] = call_data["Species"]
            self.seqwho_data[s_name]["predicted_library"] = call_data["Library"]
            self.seqwho_data[s_name]["mle"] = call_data["Maximum Likelihood Estimate"]
            self.seqwho_data[s_name]["est_read_number"] = parsed_json[s_name]["Estimated Read Number"]

            # Make MLE Matrix table for subset of predicted species and library
            mle_matrix = call_data["MLE Matrix"]
            mle_columns = [
                "human",
                "mouse",
                "amplicon",
                "bisulf",
                "wxs",
                "chip",
                "wgs",
                "dnase",
                "rnaseq",
                "atac",
                "mirnaseq",
            ]
            for i in mle_columns:
                self.seqwho_data[s_name][i] = mle_matrix[i]

            # Make a Quality Score vs count  version of the data
            self.seqwho_qualdis[s_name] = OrderedDict()
            qual_dist = parsed_json[s_name]["Quality Dist"]
            for i in range(len(qual_dist)):
                self.seqwho_qualdis[s_name][i] = qual_dist[i]

            # Make a Quality score vs position  version of the data
            self.seqwho_qualscore[s_name] = OrderedDict()
            qual_scores = parsed_json[s_name]["Quality Scores"]
            for i in range(len(qual_scores)):
                self.seqwho_qualscore[s_name][i] = qual_scores[i]

            # Make a Read Length vs count of the data
            self.seqwho_readdist[s_name] = OrderedDict()
            read_dist = parsed_json[s_name]["Read lengths"]
            for i in range(len(read_dist)):
                self.seqwho_readdist[s_name][i] = read_dist[i]

    def seqwho_general_stats_table(self):
        """Take the parsed stats from the SeqWho report and add it to the
        General Statistics table at the top of the report"""

        headers = OrderedDict()
        headers["predicted_species"] = {
            "title": "Predicted Species",
            "description": "Predicted Species",
        }
        headers["predicted_library"] = {
            "title": "Predicted Library type",
            "description": "Predicted Library type",
        }
        headers["mle"] = {
            "title": "Maximum Likelihood Estimate (Overall)",
            "description": "Overall Maximum Likelihood Estimate",
            "min": 0,
            "format": "{:,.3f}",
            "scale": False,
        }
        headers["est_read_number"] = {
            "title": "Estimated Read Number (millions)",
            "description": "M Estimated Read Number",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.general_stats_addcols(self.seqwho_data, headers)

    def seqwho_species_plot(self):
        """Take the Maximimum Likelihood Species data from the SeqWho report and make
        barchart"""

        headers = ["human", "mouse"]
        # Config for the plot
        config = {
            "id": "SeqWho_species_plot",
            "title": "SeqWho: Species",
            "ylab": "Sample",
            "ymax": 1.0,
            "stacking": None,
            "cpswitch": False,
            "cpswitch_counts_label": "Maximum Likelihood",
            "tt_decimals": 3,
        }

        return bargraph.plot(self.seqwho_data, headers, config)

    def seqwho_library_plot(self):
        """Take the Maximimum Likelihood Library data from the SeqWho report and make a
        barplot"""

        headers = ["amplicon", "bisulf", "wxs", "chip", "wgs", "dnase", "rnaseq", "atac", "mirnaseq"]
        # Config for the plot
        config = {
            "id": "SeqWho_sequencing_plot",
            "title": "SeqWho: Library",
            "ylab": "Sample",
            "ymax": 1.0,
            "stacking": None,
            "cpswitch": False,
            "cpswitch_counts_label": "Maximum Likelihood",
            "tt_decimals": 3,
        }

        return bargraph.plot(self.seqwho_data, headers, config)

    def seqwho_qualdist_plot(self):
        """Generate the Average Sequence Quality Distribution plot"""
        pconfig = {
            "id": "SeqWho_qualdis_plot",
            "title": "SeqWho: Per Sequence Quality Scores",
            "ylab": "Counts",
            "xlab": "Sequence Quality (Phred)",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>Phred {point.x}</b>: {point.y:.0f} reads",
        }

        self.add_section(
            name="Predicted Per Sequence Quality Scores",
            description="This plot shows the predicted number of reads with average quality scores.",
            plot=linegraph.plot([self.seqwho_qualdis], pconfig),
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
            description="This plot shows the predicted mean quality value across each base position in the read.",
            plot=linegraph.plot([self.seqwho_qualscore], pconfig),
        )

    def seqwho_readdist_plot(self):
        """Generate the number of reads with certain length plot"""
        pconfig = {
            "id": "seqwho_readdist_plot",
            "title": "SeqWho: Read Distribution",
            "ylab": "Counts",
            "xlab": "Length (bp)",
            "xDecimals": False,
            "ymin": 0,
            "tt_label": "<b>Base {point.x} bp</b>: {point.y} reads",
        }

        self.add_section(
            name="Predicted Read Distribution",
            description="This plot shows predicted the number of reads with certain lengths.",
            plot=linegraph.plot([self.seqwho_readdist], pconfig),
        )
