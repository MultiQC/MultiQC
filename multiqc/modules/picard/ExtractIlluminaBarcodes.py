#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard ExtractIlluminaBarcodes """

from collections import OrderedDict
import logging
import re

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard ExtractIlluminaBarcodes reports and parse their data"""

    # Set up vars
    self.picard_barcode_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/extractilluminabarcodes", filehandles=True):
        lane = None
        raw_data = []
        keys = None

        for line in f["f"]:
            # Pull lane number from input
            if "ExtractIlluminaBarcodes" in line and "LANE" in line:
                lane_search = re.findall(r"LANE(?:=|\s+)(\[?[^\s]+\]?)", line, flags=re.IGNORECASE)
                if lane_search:
                    lane = lane_search[0]
                    self.picard_barcode_metrics[lane] = {}
            if "ExtractIlluminaBarcodes" in line and "## METRICS CLASS" in line:
                keys = f["f"].readline().strip("\n").split("\t")
            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) == len(keys):
                    data = dict(zip(keys, vals))
                    data["LANE"] = lane
                    raw_data.append(data)

        for bc_data in raw_data:
            self.picard_barcode_metrics[bc_data["LANE"]][bc_data["BARCODE"]] = bc_data

    # Filter to strip out ignored sample names
    self.picard_barcode_metrics = self.ignore_samples(self.picard_barcode_metrics)

    # Stop if we didn't find anything
    if len(self.picard_barcode_metrics) == 0:
        return 0

    # Write parsed data to a file
    self.write_data_file(self.picard_barcode_metrics, "multiqc_picard_ExtractIlluminaBarcodes")

    plot_data = {}
    plot_data["per_lane"] = reads_per_lane(self.picard_barcode_metrics)
    plot_data["per_bc"] = reads_per_barcode(self.picard_barcode_metrics)

    per_lane_plot_config = {
        "id": "plot-picard-illuminabarcodemetrics-readsperlane",
        "title": "Picard ExtractIlluminaBarcodes: Reads per lane",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Matches", "ylab": "Number of Matching Reads"},
            {"name": "PF Matches", "ylab": "Number of Passing Filter Reads"},
        ],
    }
    per_barcode_plot_config = {
        "id": "plot-picard-illuminabarcodemetrics-readsperbarcode",
        "title": "Picard ExtractIlluminaBarcodes: Reads per barcode",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Matches", "ylab": "Number of Matching Reads"},
            {"name": "PF Matches", "ylab": "Number of Passing Filter Reads"},
        ],
    }

    plot_cats = [OrderedDict(), OrderedDict(), OrderedDict()]
    plot_cats[0]["READS"] = {"name": "Reads"}
    plot_cats[1]["PERFECT_MATCHES"] = {"name": "Perfect Matching Reads"}
    plot_cats[1]["ONE_MISMATCH_MATCHES"] = {"name": "One Mismatch Reads"}
    plot_cats[2]["PF_PERFECT_MATCHES"] = {"name": "Perfect Matching Passing Filter Reads"}
    plot_cats[2]["PF_ONE_MISMATCH_MATCHES"] = {"name": "Passing Filter One Mismatch Reads"}

    self.add_section(
        name="Barcode Metrics Per Lane",
        anchor="picard-illuminabarcodemetrics-perlane",
        description="""
            Indicates the number of matches (and mismatches) between the barcode reads and the actual barcodes.
            See the [Picard Documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric) for details.
        """,
        plot=bargraph.plot(
            [plot_data["per_lane"], plot_data["per_lane"], plot_data["per_lane"]],
            plot_cats,
            per_lane_plot_config,
        ),
    )

    self.add_section(
        name="Barcode Metrics Per Barcode",
        anchor="picard-illuminabarcodemetrics-perbarcode",
        description="""
            Indicates the number of matches (and mismatches) between the barcode reads and the actual barcodes.
            See the [Picard Documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric) for details.
        """,
        plot=bargraph.plot(
            [plot_data["per_bc"], plot_data["per_bc"], plot_data["per_bc"]],
            plot_cats,
            per_barcode_plot_config,
        ),
    )
    # Return the number of detected samples to the parent module
    return len(self.picard_barcode_metrics)


def reads_per_barcode(data):
    """
    Return reads/barcode for entire run
    """
    reads_per_barcode = {}
    for lane, barcodes in data.items():
        for barcode, barcode_data in barcodes.items():
            if not barcode in reads_per_barcode:
                reads_per_barcode[barcode] = {
                    "READS": 0,
                    "PERFECT_MATCHES": 0,
                    "ONE_MISMATCH_MATCHES": 0,
                    "PF_READS": 0,
                    "PF_PERFECT_MATCHES": 0,
                    "PF_ONE_MISMATCH_MATCHES": 0,
                }
            reads_per_barcode[barcode]["READS"] += int(barcode_data["READS"])
            reads_per_barcode[barcode]["PERFECT_MATCHES"] += int(barcode_data["PERFECT_MATCHES"])
            reads_per_barcode[barcode]["ONE_MISMATCH_MATCHES"] += int(barcode_data["ONE_MISMATCH_MATCHES"])
            reads_per_barcode[barcode]["PF_READS"] += int(barcode_data["PF_READS"])
            reads_per_barcode[barcode]["PF_PERFECT_MATCHES"] += int(barcode_data["PF_PERFECT_MATCHES"])
            reads_per_barcode[barcode]["PF_ONE_MISMATCH_MATCHES"] += int(barcode_data["PF_ONE_MISMATCH_MATCHES"])
    return reads_per_barcode


def reads_per_lane(data):
    """
    Return total reads/lane
    """
    reads_per_lane = {}
    for lane, barcodes in data.items():
        reads_per_lane[lane] = {
            "READS": 0,
            "PERFECT_MATCHES": 0,
            "ONE_MISMATCH_MATCHES": 0,
            "PF_READS": 0,
            "PF_PERFECT_MATCHES": 0,
            "PF_ONE_MISMATCH_MATCHES": 0,
        }
        for barcode, barcode_data in barcodes.items():
            reads_per_lane[lane]["READS"] += int(barcode_data["READS"])
            reads_per_lane[lane]["PERFECT_MATCHES"] += int(barcode_data["PERFECT_MATCHES"])
            reads_per_lane[lane]["ONE_MISMATCH_MATCHES"] += int(barcode_data["ONE_MISMATCH_MATCHES"])
            reads_per_lane[lane]["PF_READS"] += int(barcode_data["PF_READS"])
            reads_per_lane[lane]["PF_PERFECT_MATCHES"] += int(barcode_data["PF_PERFECT_MATCHES"])
            reads_per_lane[lane]["PF_ONE_MISMATCH_MATCHES"] += int(barcode_data["PF_ONE_MISMATCH_MATCHES"])

    return reads_per_lane
