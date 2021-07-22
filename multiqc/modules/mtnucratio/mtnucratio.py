#!/usr/bin/env python

""" MultiQC module to parse output from mtnucratio """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """mtnucratio module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="mtnucratio",
            anchor="mtnucratio",
            href="http://www.github.com/apeltzer/MTNucRatioCalculator",
            info="is a tool to compute mt/nuc ratios for NGS datasets.",
        )

        # Find and load any MTNUCRATIO reports
        self.mtnuc_data = dict()

        for f in self.find_log_files("mtnucratio", filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.mtnuc_data = self.ignore_samples(self.mtnuc_data)

        if len(self.mtnuc_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.mtnuc_data)))

        # Write parsed report data to a file
        self.write_data_file(self.mtnuc_data, "multiqc_mtnucratio")

        # Basic Stats Table
        self.mtnucratio_general_stats_table()

    # Parse our nice little JSON file
    def parseJSON(self, f):
        """Parse the JSON output from mtnucratio and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
            if "metrics" not in parsed_json and "metadata" not in parsed_json:
                log.warning("No MTNUCRATIO JSON: '{}'".format(f["fn"]))
                return None
        except json.JSONDecodeError as e:
            log.warning("Could not parse mtnucratio JSON: '{}'".format(f["fn"]))
            log.debug(e)
            return None

        # Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json["metadata"]["sample_name"], f)
        self.add_data_source(f, s_name)

        metrics_dict = parsed_json["metrics"]

        # Add all in the main data_table
        self.mtnuc_data[s_name] = metrics_dict

    def mtnucratio_general_stats_table(self):
        """Take the parsed stats from the mtnucratio report and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["mt_cov_avg"] = {
            "title": "MT genome coverage",
            "description": "Average coverage (X) on mitochondrial genome.",
            "min": 0,
            "scale": "OrRd",
            "suffix": " X",
            "hidden": True,
        }
        headers["nuc_cov_avg"] = {
            "title": "Genome coverage",
            "description": "Average coverage (X) on nuclear genome.",
            "min": 0,
            "scale": "GnBu",
            "suffix": " X",
            "hidden": True,
        }
        headers["mt_nuc_ratio"] = {
            "title": "MT to Nuclear Ratio",
            "description": "Mitochondrial to nuclear reads ratio (MTNUC)",
            "min": 0,
            "max": 100,
            "scale": "RdYlGrn-rev",
        }
        headers["nucreads"] = {
            "title": "{} Genome reads".format(config.read_count_prefix),
            "description": "Reads on the nuclear genome ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "scale": "BuPu",
            "hidden": True,
        }
        headers["mtreads"] = {
            "title": "{} MT genome reads".format(config.read_count_prefix),
            "description": "Reads on the mitochondrial genome ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "scale": "OrRd",
            "hidden": True,
        }

        self.general_stats_addcols(self.mtnuc_data, headers)
