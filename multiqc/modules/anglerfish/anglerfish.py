#!/usr/bin/env python

""" MultiQC module to parse output from Anglerfish """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, linegraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Anglerfish module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="anglerfish",
            anchor="anglerfish",
            href="https://github.com/remiolsen/anglerfish",
            info="A tool to assess Illumina libraries sequenced on Oxford Nanopore for the purpose of quality control.",
        )

        # Find and load any anglerfish reports
        self.anglerfish_data = dict()

        for f in self.find_log_files("anglerfish", filehandles=True):
            self.parse_anglerfish_json(f)

        # Filter to strip out ignored sample names
        self.anglerfish_data = self.ignore_samples(self.anglerfish_data)

        # Stop execution of the data if no anglerfish data is found.
        if len(self.anglerfish_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.anglerfish_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.anglerfish_data, "multiqc_anglerfish")

        # General Stats Table
        self.anglerfish_general_stats_table()

        # Paf Statistics bar plot
        self.add_section(
            name="Paf Statistics",
            anchor="anglerfish-paf-statistics",
            description="Paf Statistics of sampled reads.",
            plot=self.anglerfish_paf_stats_chart(),
        )
        # Sample Stats -mean scatter plot
        self.add_section(
            name="Sample Statistics",
            anchor="anglerfish-sample-statistics",
            description="Outliers for the sample statistics",
            plot=self.anglerfish_sample_stats_chart(),
        )
        # Undetermined plot
        # TODO: error handle if no undetermined exists?
        self.add_section(
            name="Undetermined indexes",
            anchor="anglerfish-undetermined-indexes",
            description="Showcases undetermined indexes",
            plot=self.anglerfish_undetermined_index_chart(),
        )

    def parse_anglerfish_json(self, f):
        """Parse the JSON output from anglerfish and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except:
            log.warning("Could not parse anglerfish JSON: '{}'".format(f["fn"]))
            return None

        # Fetch a sample name from the command
        s_name = f["s_name"]
        self.add_data_source(f, s_name)
        self.anglerfish_data[s_name] = {}

        # Parse paf statistics
        index = 0
        for k in parsed_json["paf_stats"]:
            key_list = k.keys()
            for key in key_list:
                try:
                    self.anglerfish_data[s_name][key + "_{}".format(index)] = float(k[key][0])
                except KeyError:
                    log.debug("'" + key + "' key missing in anglerfish json '{}'".format(f["fn"]))
            index += 1
        # Save total index amount, for plotting purpose
        self.anglerfish_data[s_name]["paf_stats_amount"] = index

        # Parse Sample Reads
        index = 0
        for k in parsed_json["sample_stats"]:
            key_list = k.keys()
            for key in key_list:
                if key != "sample_name":
                    try:
                        self.anglerfish_data[s_name][key + "_{}".format(index)] = float(k[key])
                    except:
                        log.debug("'" + key + "' missing in anglerfish json: '{}'".format(f["fn"]))
                else:
                    try:
                        self.anglerfish_data[s_name]["sample_name_{}".format(index)] = k["sample_name"]
                    except:
                        log.debug("'sample_name' key missing in anglerfish json: '{}'".format(f["fn"]))
            index += 1
        self.anglerfish_data[s_name]["sample_stats_amount"] = index

        # Parse Undetermined
        index = 0
        for k in parsed_json["undetermined"]:
            self.anglerfish_data[s_name]["undetermined_count_{}".format(index)] = float(k["count"])
            self.anglerfish_data[s_name]["undetermined_index_{}".format(index)] = k["undetermined_index"]
            index += 1
        self.anglerfish_data[s_name]["undetermined_amount"] = index

    # TODO: General stats table
    def anglerfish_general_stats_table(self):
        """Add Anglerfish statistics to the general statistics table"""
        headers = OrderedDict()
        headers["librarys"] = {
            "title": "Librarys",
            "description": "Amount of librarys used",
            # TODO: Rest of info
        }
        # TODO: Add headers
        # library header?

        self.general_stats_addcols(self.anglerfish_data, headers, "anglerfish")

    def anglerfish_paf_stats_chart(self):
        """Generate the Anglerfish Paf stats plot"""

        # Data structure for the Paf stat plot
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["paf_stats_amount"]
            for i in range(index):
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)] = {}
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)][
                    "aligned reads matching both I7 and I5 adaptor"
                ] = self.anglerfish_data[s_name]["aligned reads matching both I7 and I5 adaptor_{}".format(i)]
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)][
                    "aligned reads matching multiple I7/I5 adaptor pairs"
                ] = self.anglerfish_data[s_name]["aligned reads matching multiple I7/I5 adaptor pairs_{}".format(i)]
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)][
                    "aligned reads matching only I7 or I5 adaptor"
                ] = self.anglerfish_data[s_name]["aligned reads matching only I7 or I5 adaptor_{}".format(i)]
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)][
                    "aligned reads with uncategorized alignments"
                ] = self.anglerfish_data[s_name]["aligned reads with uncategorized alignments_{}".format(i)]
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)]["input_reads"] = self.anglerfish_data[s_name][
                    "input_reads_{}".format(i)
                ]
                data["{s} Paf Stats, {i}".format(s=s_name, i=i)][
                    "reads aligning to adaptor sequences"
                ] = self.anglerfish_data[s_name]["reads aligning to adaptor sequences_{}".format(i)]
        config = {
            "id": "Anglerfish_paf_plot",
            "cpswitch": False,
            "title": "Anglerfish: Paf Plot",
            "stacking": None,
        }
        # return bargraph.plot(data, keys, config)
        return bargraph.plot(data, None, config)

    def anglerfish_sample_stats_chart(self):
        """Generate Sample Stats Scatter Plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["sample_stats_amount"]

            for i in range(index):
                sample_name = self.anglerfish_data[s_name]["sample_name_{}".format(i)]
                data["std_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)] = {}
                data["mean_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)] = {}

                data["std_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["x"] = self.anglerfish_data[
                    s_name
                ]["#reads_{}".format(i)]
                data["mean_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["x"] = self.anglerfish_data[
                    s_name
                ]["#reads_{}".format(i)]

                data["std_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["y"] = self.anglerfish_data[
                    s_name
                ]["std_read_len_{}".format(i)]
                data["mean_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["y"] = self.anglerfish_data[
                    s_name
                ]["mean_read_len_{}".format(i)]

                data["std_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["color"] = "#1f78b4"
                data["mean_read_len Sample: {i}, {s}".format(i=sample_name, s=s_name)]["color"] = "#33a02c"
            config = {
                "id": "sample_stats_scatter_plot",
                "title": "Anglerfish: Sample Stats",
                "xlab": "#reads",
                "ylab": "read_len",
                "ymin": 0,
                "xmin": 0,
            }
        return scatter.plot(data, config)

    def anglerfish_undetermined_index_chart(self):
        """Generate Undetermined indexes Bar Plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["undetermined_amount"]
            for i in range(index):
                undetermined_index = self.anglerfish_data[s_name]["undetermined_index_{}".format(i)]
                data["{s}: {u_i}".format(s=s_name, u_i=undetermined_index)] = {}
                data["{s}: {u_i}".format(s=s_name, u_i=undetermined_index)][undetermined_index] = self.anglerfish_data[
                    s_name
                ]["undetermined_count_{}".format(i)]
        config = {
            "id": "Anglerfish_undetermined_index_plot",
            "cpswitch": False,
            "title": "Anglerfish: Undetermined Indexes",
        }
        return bargraph.plot(data, None, config)
