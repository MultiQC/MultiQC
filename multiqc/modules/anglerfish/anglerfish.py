#!/usr/bin/env python

""" MultiQC module to parse output from Anglerfish """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, scatter, beeswarm, table
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
            name="Anglerfish",
            anchor="Anglerfish",
            href="https://github.com/remiolsen/anglerfish",
            info="A tool to assess Illumina libraries sequenced on Oxford Nanopore for the purpose of quality control.",
        )

        # Find and load any anglerfish reports
        self.anglerfish_data = dict()
        self.anglerfish_gst = dict()

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

        # Sample Stats Read length table/beeswarm plot
        self.add_section(
            name="Read Lengths Summary",
            anchor="Anglerfish-sample-statistics",
            description="The mean read length and the standard deviation for each sample.",
            plot=self.anglerfish_sample_stats(),
        )
        # Undetermined indexes plot
        self.add_section(
            name="Undetermined indexes",
            anchor="Anglerfish-undetermined-indexes",
            description="",
            plot=self.anglerfish_undetermined_index_chart(),
        )

    def parse_anglerfish_json(self, f):
        """Parse the JSON output from Anglerfish and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except:
            log.warning("Could not parse Anglerfish JSON: '{}'".format(f["fn"]))
            return None

        # Fetch a sample name from the command
        s_name = f["s_name"]
        self.add_data_source(f, s_name)
        self.anglerfish_data[s_name] = {}

        # Parse Sample Stats
        index = 0
        for k in parsed_json["sample_stats"]:
            key_list = k.keys()
            for key in key_list:
                if key != "sample_name":
                    self.anglerfish_data[s_name][key + "_{}".format(index)] = float(k[key])
                else:
                    self.anglerfish_data[s_name]["sample_name_{}".format(index)] = k["sample_name"]
            index += 1
        self.anglerfish_data[s_name]["sample_stats_amount"] = index

        # Parse Undetermined Indexes
        index = 0
        total_count = 0
        for k in parsed_json["undetermined"]:
            self.anglerfish_data[s_name]["undetermined_count_{}".format(index)] = float(k["count"])
            total_count += float(k["count"])
            self.anglerfish_data[s_name]["undetermined_index_{}".format(index)] = k["undetermined_index"]
            index += 1
        self.anglerfish_data[s_name]["undetermined_amount"] = index

        # Parse for general stat table
        ## Multiple sample names per file requires dict where thr first key is not file name
        total_read = parsed_json["paf_stats"][0]["input_reads"][0]
        for k in parsed_json["sample_stats"]:
            key = k["sample_name"]
            reads = float(k["#reads"])
            self.anglerfish_gst[key] = {}
            self.anglerfish_gst["undetermined"] = {}
            self.anglerfish_gst[key]["#reads"] = reads
            self.anglerfish_gst[key]["mean_read_len"] = float(k["mean_read_len"])
            self.anglerfish_gst[key]["std_read_len"] = float(k["std_read_len"])
            try:
                self.anglerfish_gst[key]["library"] = float((reads / total_read) * 100)
                self.anglerfish_gst["undetermined"]["library"] = float((total_count / total_read) * 100)
            except (ZeroDivisionError):
                log.debug("zero input reads - library = 0")
                self.anglerfish_gst[key]["library"] = float(0)
                self.anglerfish_gst["undetermined"]["library"] = float(0)

    # General stats table
    def anglerfish_general_stats_table(self):
        """Add Anglerfish statistics to the general statistics table"""
        headers = OrderedDict()
        headers["library"] = {
            "namespace": "Anglerfish",
            "title": "library",
            "description": "Amount of library used. Calculated from a samples #reads divided by total number of input reads",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn-rev",
            "suffix": " %",
        }

        headers["#reads"] = {
            "namespace": "Anglerfish",
            "title": "#reads",
            "description": "Amount of reads",
            "min": 0,
            "scale": "RdYlGn-rev",
        }
        headers["mean_read_len"] = {
            "namespace": "Anglerfish",
            "title": "mean_read_len",
            "description": "Mean read length",
            "min": 0,
            "scale": "RdYlGn-rev",
        }
        headers["std_read_len"] = {
            "namespace": "Anglerfish",
            "title": "std_read_len",
            "description": "Standard deviation read length",
            "min": 0,
            "scale": "RdYlGn-rev",
        }

        self.general_stats_addcols(self.anglerfish_gst, headers, "anglerfish")

    def anglerfish_sample_stats(self):
        """Generate plot for read length from sample stats.
        For >10 samples: generate table plot
        for >= 10 samples: generate beeswarm plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["sample_stats_amount"]

            for i in range(index):
                sample_name = self.anglerfish_data[s_name]["sample_name_{}".format(i)]
                data["Sample: {}".format(sample_name)] = {}
                data["Sample: {}".format(sample_name)]["Mean"] = self.anglerfish_data[s_name][
                    "mean_read_len_{}".format(i)
                ]
                data["Sample: {}".format(sample_name)]["Standard Deviation"] = self.anglerfish_data[s_name][
                    "std_read_len_{}".format(i)
                ]

        # Plot table if less than 10 samples exist, beeswarm if more
        if index < 10:
            return table.plot(data)
        return beeswarm.plot(data)

    def anglerfish_undetermined_index_chart(self):
        """Generate Undetermined indexes Bar Plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["undetermined_amount"]
            for i in range(index):
                undetermined_index = self.anglerfish_data[s_name]["undetermined_index_{}".format(i)]
                data["{}".format(undetermined_index)] = {}
                data["{}".format(undetermined_index)][undetermined_index] = self.anglerfish_data[s_name][
                    "undetermined_count_{}".format(i)
                ]
        config = {
            "id": "Anglerfish_undetermined_index_plot",
            "cpswitch": False,
            "title": "Anglerfish: Undetermined Indexes",
            "tt_percentages": False,
        }
        return bargraph.plot(data, None, config)
