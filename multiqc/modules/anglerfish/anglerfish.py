#!/usr/bin/env python

""" MultiQC module to parse output from Anglerfish """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, scatter
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
                    self.anglerfish_data[s_name][key + "_P_{}".format(index)] = float((k[key][1]) * 100)
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
        total_count = 0
        for k in parsed_json["undetermined"]:
            self.anglerfish_data[s_name]["undetermined_count_{}".format(index)] = float(k["count"])
            total_count += float(k["count"])
            self.anglerfish_data[s_name]["undetermined_index_{}".format(index)] = k["undetermined_index"]
            index += 1
        self.anglerfish_data[s_name]["undetermined_amount"] = index

        # Parse Library
        total_read = parsed_json["paf_stats"][0]["input_reads"][0]
        for k in parsed_json["sample_stats"]:
            key = k["sample_name"]
            reads = float(k["#reads"])
            self.anglerfish_gst[key] = {}
            self.anglerfish_gst["undetermined"] = {}
            self.anglerfish_gst[key]["library"] = float((reads / total_read) * 100)
            self.anglerfish_gst[key]["#reads"] = reads
        self.anglerfish_gst["undetermined"]["library"] = float((total_count / total_read) * 100)

    # General stats table
    def anglerfish_general_stats_table(self):
        """Add Anglerfish statistics to the general statistics table"""
        headers = OrderedDict()
        headers["library"] = {
            "title": "library",
            "description": "Amount of library used",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn-rev",
            "suffix": " %",
        }

        headers["#reads"] = {
            "title": "#reads",
            "description": "Amount of reads",
            "min": 0,
            "scale": "RdYlGn-rev",
        }

        self.general_stats_addcols(self.anglerfish_gst, headers, "anglerfish")

    def anglerfish_paf_stats_chart(self):
        """Generate the Anglerfish Paf stats plot"""

        # Data structure for the Paf stat plot
        ## Data Structure for grouped amounts
        dataG = {}
        ## Data structure for single amounts
        dataS = {}
        ## Data Structure for grouped percentages
        dataG_P = {}
        ## Data structure for single percentages
        dataS_P = {}
        # Keys
        key_list = [
            "aligned reads matching both I7 and I5 adaptor",
            "aligned reads matching multiple I7/I5 adaptor pairs",
            "aligned reads matching only I7 or I5 adaptor",
            "aligned reads with uncategorized alignments",
            "input_reads",
            "reads aligning to adaptor sequences",
        ]
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["paf_stats_amount"]
            if index == 1:
                self.group_data(dataG, s_name, key_list, None, None)
                self.group_data(dataG_P, s_name, key_list, True, None)
                self.single_data(dataS, s_name, key_list, None, None)
                self.single_data(dataS_P, s_name, key_list, True, None)
            else:
                for i in range(index):
                    self.group_data(dataG, s_name, key_list, None, i)
                    self.group_data(dataG_P, s_name, key_list, True, i)
                    self.single_data(dataS, s_name, key_list, None, i)
                    self.single_data(dataS_P, s_name, key_list, True, i)

        config = {
            "id": "Anglerfish_paf_plot",
            "cpswitch": False,
            "data_labels": [
                {"name": "Amount in Group", "ylab": "Amount matched"},
                {"name": "Amount Single", "ylab": "Amount matched"},
                {"name": "Percentages in Group", "ylab": "Percent", "ymax": 100},
                {"name": "Percentages Single", "ylab": "Percent", "ymax": 100},
            ],
            "title": "Anglerfish: Paf Plot",
            "stacking": None,
            "tt_decimals": 2,
            "tt_percentages": False,
        }
        return bargraph.plot([dataG, dataS, dataG_P, dataS_P], None, config)

    def group_data(self, data, s_name, key_list, percent=None, index=None):
        """Make a data structure with all keys in a group"""
        keyN = "Paf stats"
        suffix = "_0"
        if index != None:
            keyN = "Paf Stats, {}".format(index)
            suffix = "_{}".format(index)
        data[keyN] = {}
        if percent != None:
            suffix = "_P" + suffix
        for key in key_list:
            data[keyN][key] = self.anglerfish_data[s_name][key + suffix]

    def single_data(self, data, s_name, key_list, percent=None, index=None):
        """Make data structure with all keys as own categories"""
        suffix = "_0"
        if index != None:
            suffix = "_{}".format(index)
        if percent != None:
            suffix = "_P" + suffix
        for key in key_list:
            keyN = key
            if index != None:
                keyN = key + "_{}".format(index)
            data[keyN] = {}
            data[keyN][key] = self.anglerfish_data[s_name][key + suffix]

    def anglerfish_sample_stats_chart(self):
        """Generate Sample Stats Scatter Plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["sample_stats_amount"]

            for i in range(index):
                sample_name = self.anglerfish_data[s_name]["sample_name_{}".format(i)]
                data["Sample: {}".format(sample_name)] = {}
                data["Sample: {}".format(sample_name)]["x"] = self.anglerfish_data[s_name]["std_read_len_{}".format(i)]
                data["Sample: {}".format(sample_name)]["y"] = self.anglerfish_data[s_name]["mean_read_len_{}".format(i)]
            config = {
                "id": "sample_stats_scatter_plot",
                "title": "Anglerfish: Sample Stats",
                "xlab": "std_read_len",
                "ylab": "mean_read_len",
            }
        return scatter.plot(data, config)

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
