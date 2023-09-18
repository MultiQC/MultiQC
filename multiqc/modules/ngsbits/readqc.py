#!/usr/bin/env python

""" MultiQC module to parse output from ngs-bits ReadQC tool """

import logging
import xml.etree.cElementTree

from multiqc import config
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


def check_paired_end(qcml_contents):
    """Check if both R1 and R2 are present in input files."""
    found_r1 = False
    found_r2 = False
    root = xml.etree.cElementTree.fromstring(qcml_contents)

    for el in root.findall(".//{http://www.prime-xs.eu/ms/qcml}metaDataParameter"):
        if el.attrib["name"] == "source file":
            if "R1" in el.attrib["value"]:
                found_r1 = True
            if "R2" in el.attrib["value"]:
                found_r2 = True

    return found_r1 and found_r2


def parse_reports(self):
    """Find ngs-bits ReadQC reports and parse their data"""

    # Set up vars
    self.readqc = dict()
    self.readqc_keys = dict()

    for f in self.find_log_files("ngsbits/readqc"):
        d = self.parse_qcml_by(f["f"], "qualityParameter")
        is_pe = check_paired_end(f["f"])

        if len(d[0]) > 0:
            if f["s_name"] in self.readqc:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
            self.add_data_source(f, section="readqc")
            self.readqc[f["s_name"]] = d[0]
            self.readqc_keys.update(d[1])
            self.readqc[f["s_name"]].update(
                {
                    "cluster count": self.readqc[f["s_name"]]["read count"] / (1 + is_pe),
                    "paired-end": "yes" if is_pe else "no",
                }
            )

    # Filter to strip out ignored sample names
    self.readqc = self.ignore_samples(self.readqc)

    if len(self.readqc) > 0:
        # Convert numbers given in megabases to bases
        self.readqc_keys["bases sequenced"] = ("Bases sequenced in total.", "")
        for _, kv in self.readqc.items():
            kv["bases sequenced"] = kv["bases sequenced (MB)"] * 1e6
            kv.pop("bases sequenced (MB)")

        # Add cluster count
        self.readqc_keys["cluster count"] = ("Clusters sequenced in total.", "")
        self.readqc_keys["paired-end"] = ("Whether input files were paired-end sequences.", "")

        # Improve table headers
        self.readqc_keys_table = {
            key: {"title": key, "description": value[0]} for key, value in self.readqc_keys.items()
        }

        self.readqc_keys_table["read count"].update(
            {
                "suffix": config.read_count_prefix,
                "format": "{:,.2f}",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "Purples",
                "placement": 10,
            }
        )
        self.readqc_keys_table["cluster count"].update(
            {
                "suffix": config.read_count_prefix,
                "format": "{:,.2f}",
                "modify": lambda x: x * config.read_count_multiplier,
                "scale": "Purples",
                "placement": 20,
            }
        )
        self.readqc_keys_table["bases sequenced"].update(
            {
                "suffix": config.base_count_prefix,
                "format": "{:,.2f}",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
                "scale": "Blues",
                "placement": 30,
            }
        )
        self.readqc_keys_table["gc content %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Spectral", "placement": 40}
        )
        self.readqc_keys_table["Q20 read %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Reds", "placement": 50}
        )
        self.readqc_keys_table["Q30 base %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Oranges", "placement": 60}
        )
        self.readqc_keys_table["read length"].update(
            {"suffix": "bp", "format": "{:,.0f}", "scale": "Greens", "placement": 70}
        )
        self.readqc_keys_table["no base call %"].update(
            {"suffix": "%", "format": "{:,.2f}", "floor": 1, "scale": "BuGn", "placement": 80}
        )

        # Write to file
        self.write_data_file(self.readqc, "multiqc_ngsbits_readqc")

        # overview table with all values
        self.add_section(
            name="ReadQC",
            anchor="ngsbits-readqc",
            description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md" target="_blank">ReadQC</a>'
            " calculates QC metrics on unprocessed NGS reads.",
            plot=table.plot(self.readqc, self.readqc_keys_table, pconfig={"namespace": "ngsbits_readqc"}),
        )

    return len(self.readqc)
