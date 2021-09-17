#!/usr/bin/env python

""" MultiQC module to parse output from ngs-bits MappingQC tool """

import logging

from multiqc import config
from multiqc.plots import table

import xml.etree.cElementTree

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find ngs-bits MappingQC reports and parse their data"""

    # Set up vars
    self.mappingqc = dict()
    self.mappingqc_keys = dict()

    for f in self.find_log_files("ngsbits/mappingqc"):
        d = self.parse_qcml_by(f["f"], "qualityParameter")

        if len(d[0]) > 0:
            if f["s_name"] in self.mappingqc:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
            self.add_data_source(f, section="mappingqc")
            self.mappingqc[f["s_name"]] = d[0]
            self.mappingqc_keys.update(d[1])

    # Filter to strip out ignored sample names
    self.mappingqc = self.ignore_samples(self.mappingqc)

    if len(self.mappingqc) > 0:

        # Convert numbers given in megabases to bases
        self.mappingqc_keys["bases usable"] = ("Bases usable in total.", "")
        for _, kv in self.mappingqc.items():
            kv["bases usable"] = kv["bases usable (MB)"] * 1e6
            kv.pop("bases usable (MB)")

        # Improve table headers
        self.mappingqc_keys_table = {
            key: {"title": key, "description": value[0]} for key, value in self.mappingqc_keys.items()
        }

        self.mappingqc_keys_table["trimmed base %"].update(
            {"suffix": "%", "format": "{:,.2f}", "floor": 1, "scale": "PuBu"}
        )
        self.mappingqc_keys_table["clipped base %"].update(
            {"suffix": "%", "format": "{:,.2f}", "floor": 1, "scale": "PuRd"}
        )
        self.mappingqc_keys_table["mapped read %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Reds"}
        )
        self.mappingqc_keys_table["bases usable"].update(
            {
                "suffix": config.base_count_prefix,
                "format": "{:,.2f}",
                "modify": lambda x: x * config.base_count_multiplier,
                "scale": "Greens",
            }
        )
        # always available, even without target file
        self.mappingqc_keys_table["on-target read %"].update(
            {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "Purples"}
        )

        # only available if duplicates marked
        try:
            self.mappingqc_keys_table["duplicate read %"].update(
                {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "YlOrRd"}
            )
        except KeyError:
            pass

        # only available if paired-end
        try:
            self.mappingqc_keys_table["properly-paired read %"].update(
                {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "GnBu"}
            )
            self.mappingqc_keys_table["insert size"].update({"suffix": "bp", "format": "{:,.2f}", "scale": "RdYlGn"})
        except KeyError:
            pass

        # only available if human
        try:
            self.mappingqc_keys_table["SNV allele frequency deviation"].update(
                {"suffix": "", "format": "{:,.2f}", "floor": 0, "ceiling": 10, "minRange": 10, "scale": "Greys"}
            )
        except KeyError:
            pass

        # only available if target file provided
        coverage_values = (10, 20, 30, 50, 100, 200, 500)
        try:
            self.mappingqc_keys_table["target region read depth"].update({"suffix": "x", "format": "{:,.2f}"})
            for x in coverage_values:
                self.mappingqc_keys_table["target region {:d}x %".format(x)].update(
                    {"suffix": "%", "format": "{:,.2f}", "max": 100, "scale": "YlGn"}
                )
        except KeyError:
            pass

        # Write to file
        self.write_data_file(self.mappingqc, "multiqc_ngsbits_mappingqc")

        # overview table with all values
        self.add_section(
            name="MappingQC",
            anchor="ngsbits-mappingqc",
            description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md" target="_blank">MappingQC</a>'
            " calculates QC metrics on mapped NGS reads.",
            plot=table.plot(self.mappingqc, self.mappingqc_keys_table, pconfig={"namespace": "ngsbits_mappingqc"}),
        )

    return len(self.mappingqc)
