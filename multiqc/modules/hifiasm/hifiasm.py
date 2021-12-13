#!/usr/bin/env python

""" MultiQC module to parse output from HiFiasm """

import json
import logging
import re

from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="HiFiasm",
            anchor="hifiasm",
            href="https://github.com/chhylp123/hifiasm",
            info=": a haplotype-resolved assembler for accurate Hifi reads",
            doi="10.1038/s41592-020-01056-5",
        )

        # To store the mod data
        self.hifiasm_data = dict()
        self.parse_hifiasm_log_files()
        self.hifiasm_data = self.ignore_samples(self.hifiasm_data)

        # If we found no data
        if not self.hifiasm_data:
            raise UserWarning
        log.info("Found {} reports".format(len(self.hifiasm_data)))

        self.write_data_file(self.hifiasm_data, "multiqc_hifiasm_report")
        self.add_sections()

    def parse_hifiasm_log_files(self):
        for f in self.find_log_files("hifiasm", filehandles=True):
            filename = f["s_name"]
            if filename in self.hifiasm_data:
                log.debug(f"Duplicate sample name found! Overwriting: {filename}")
            data = self.extract_kmer_graph(f["f"])
            self.hifiasm_data[filename] = data

    def add_sections(self):
        # Plot configuration
        config = {
            "id": "hifiasm-kmr-graph",
            "title": "HiFiasm: kmer graph",
            "ylab": "Count of kmer occurrence",
            "xlab": "Kmer occurrence",
            "logswitch": True,
            "logswitch_active": True,
        }

        self.add_section(
            name="HiFiasm kmer graph",
            anchor="hifiasm-kmer-graph",
            description="Kmer counts in the input data",
            plot=linegraph.plot(self.hifiasm_data, config),
        )

    def extract_kmer_graph(self, fin):
        """Extract the kmer graph from file in"""
        data = dict()

        found_histogram = False

        for line in fin:
            if line.startswith("[M::ha_hist_line]"):
                found_histogram = True
                spline = line.strip().split()
                # Occurrence of kmer
                occurrence = spline[1][:-1]
                # Special case
                if occurrence == "rest":
                    continue
                # Count of the occurrence
                count = int(spline[3])
                data[int(occurrence)] = count
            # If we are no longer in the histogram
            elif found_histogram:
                return data
