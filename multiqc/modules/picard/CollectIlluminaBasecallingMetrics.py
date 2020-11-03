#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard CollectIlluminaBasecallingMetrics """

from collections import OrderedDict
import logging
import math
import os
import re

from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard CollectIlluminaBasecallingMetrics reports and parse their data """

    # Set up vars
    self.picard_basecalling_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(
        "picard/collectilluminabasecallingmetrics", filehandles=True
    ):
        raw_data = []
        keys = None
        for line in f["f"]:
            if "IlluminaBasecallingMetrics" and "## METRICS CLASS" in line:
                keys = f["f"].readline().strip("\n").split("\t")
            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) == len(keys):
                    raw_data.append(dict(zip(keys, vals)))
        parsed_data = {}
        for metrics_line in raw_data:

            lane = metrics_line["LANE"]
            if lane not in self.picard_basecalling_metrics:
                self.picard_basecalling_metrics[lane] = {}

            if not metrics_line["MOLECULAR_BARCODE_SEQUENCE_1"]:
                metrics_line["MOLECULAR_BARCODE_SEQUENCE_1"] = "TOTAL"
            barcode = metrics_line["MOLECULAR_BARCODE_SEQUENCE_1"]
            if barcode not in self.picard_basecalling_metrics[lane]:
                self.picard_basecalling_metrics[lane][barcode] = {}
            self.picard_basecalling_metrics[lane][barcode] = metrics_line

    if len(self.picard_basecalling_metrics) > 0:

        # Write parsed data to a file
        self.write_data_file(
            self.picard_basecalling_metrics, "multiqc_picard_IlluminaBasecallingMetrics"
        )

    # Return the number of detected samples to the parent module
    return len(self.picard_basecalling_metrics)
