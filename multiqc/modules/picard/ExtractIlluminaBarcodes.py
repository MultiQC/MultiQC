""" MultiQC submodule to parse output from Picard ExtractIlluminaBarcodes """

import logging
from collections import defaultdict

from multiqc.modules.picard import util
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard ExtractIlluminaBarcodes reports and parse their data"""

    data_by_lane = defaultdict(dict)

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/extractilluminabarcodes", filehandles=True):
        # Sample name from input file name by default
        lane = f["s_name"]
        keys = None

        for line in f["f"]:
            maybe_lane_name = util.extract_sample_name(
                module, line, f, picard_tool="ExtractIlluminaBarcodes", picard_opt="LANE"
            )
            if maybe_lane_name:
                # Starts information for a new sample
                lane = maybe_lane_name
                keys = None

            if util.is_line_right_before_table(line, picard_class=["ExtractIlluminaBarcodes", "BarcodeMetric"]):
                keys = f["f"].readline().strip("\n").split("\t")
                module.add_data_source(f, s_name=lane, section="ExtractIlluminaBarcodes")

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys):
                    keys = None
                    continue

                data = dict(zip(keys, vals))
                data["LANE"] = lane
                data_by_lane[lane][data["BARCODE"]] = data

    data_by_lane = module.ignore_samples(data_by_lane)
    if len(data_by_lane) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_lane, f"multiqc_{module.anchor}_ExtractIlluminaBarcodes")

    plot_data = {"per_lane": reads_per_lane(data_by_lane), "per_bc": reads_per_barcode(data_by_lane)}

    per_lane_plot_config = {
        "id": f"plot-{module.anchor}-illuminabarcodemetrics-readsperlane",
        "title": f"{module.name} ExtractIlluminaBarcodes: Reads per lane",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Matches", "ylab": "Number of Matching Reads"},
            {"name": "PF Matches", "ylab": "Number of Passing Filter Reads"},
        ],
    }
    per_barcode_plot_config = {
        "id": f"plot-{module.anchor}-illuminabarcodemetrics-readsperbarcode",
        "title": f"{module.name} ExtractIlluminaBarcodes: Reads per barcode",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Matches", "ylab": "Number of Matching Reads"},
            {"name": "PF Matches", "ylab": "Number of Passing Filter Reads"},
        ],
    }

    plot_cats = [dict(), dict(), dict()]
    plot_cats[0]["READS"] = {"name": "Reads"}
    plot_cats[1]["PERFECT_MATCHES"] = {"name": "Perfect Matching Reads"}
    plot_cats[1]["ONE_MISMATCH_MATCHES"] = {"name": "One Mismatch Reads"}
    plot_cats[2]["PF_PERFECT_MATCHES"] = {"name": "Perfect Matching Passing Filter Reads"}
    plot_cats[2]["PF_ONE_MISMATCH_MATCHES"] = {"name": "Passing Filter One Mismatch Reads"}

    module.add_section(
        name="Barcode Metrics Per Lane",
        anchor=f"{module.anchor}-illuminabarcodemetrics-perlane",
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

    module.add_section(
        name="Barcode Metrics Per Barcode",
        anchor=f"{module.anchor}-illuminabarcodemetrics-perbarcode",
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
    return len(data_by_lane)


def reads_per_barcode(data):
    """
    Return reads/barcode for entire run
    """
    reads_per_barcode = {}
    for lane, barcodes in data.items():
        for barcode, barcode_data in barcodes.items():
            if barcode not in reads_per_barcode:
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
