#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard CollectIlluminaLaneMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


def lane_metrics_table(data):
    headers = OrderedDict()
    headers["CLUSTER_DENSITY"] = {
        "title": "Cluster Density",
        "description": "The number of clusters per unit area on the this lane (cluster / mm^2`)",
        "scale": "Greens",
    }
    headers["TYPE_NAME"] = {
        "title": "Read Number",
        "description": "Defines an Illumina template read number (first or second)",
        "modify": lambda x: x.lower(),
    }
    headers["PREPHASING_APPLIED"] = {
        "title": "Prephasing Applied",
        "description": "Median pre-phasing value across all tiles in a lane, applied to the first and second template reads",
        "scale": "BuPu",
        "max": 1,
    }
    headers["PHASING_APPLIED"] = {
        "title": "Phasing Applied",
        "description": "Median phasing value across all tiles in a lane, applied to the first and second template reads",
        "scale": "BuPu",
        "max": 1,
    }

    table_config = {
        "id": "picard-illumina-lane-metrics-table",
        "namespace": "Picard",
        "table_title": "Picard Illumina Lane Metrics",
    }
    tdata = {}
    for run_name, run in data.items():
        for lane_number, lane in run.items():
            tdata[f"{run_name} - L{lane_number}"] = lane
    return table.plot(tdata, headers, table_config)


def parse_reports(self):
    """Find Picard CollectIlluminaLaneMetrics reports and parse their data"""

    # Set up vars
    self.picard_lane_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/collectilluminalanemetrics", filehandles=True):
        run_name = None
        parsed_data = []
        keys = []

        for line in f["f"]:
            if "CollectIlluminaLaneMetrics" in line and "RUN_DIRECTORY" in line:
                run_dir_search = re.search(r"RUN_DIRECTORY(?:=|\s+)(\[?[^\s]+\]?)", line, flags=re.IGNORECASE)
                run_name = os.path.basename(run_dir_search[0])
                if not run_name in self.picard_lane_metrics:
                    self.picard_lane_metrics[run_name] = {}

            if run_name is not None:
                if ("IlluminaLaneMetrics" in line or "IlluminaPhasingMetrics" in line) and "## METRICS CLASS" in line:
                    keys = f["f"].readline().strip("\n").split("\t")
                elif keys:
                    vals = line.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        parsed_data.append(dict(zip(keys, vals)))

                for d in parsed_data:
                    lane = d["LANE"]
                    if not lane in self.picard_lane_metrics[run_name]:
                        self.picard_lane_metrics[run_name][lane] = {}
                    self.picard_lane_metrics[run_name][lane].update(d)

    # Filter to strip out ignored sample names
    self.picard_lane_metrics = self.ignore_samples(self.picard_lane_metrics)

    if len(self.picard_lane_metrics) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_lane_metrics, "multiqc_picard_IlluminaLaneMetrics")

        self.add_section(
            name="Lane Metrics",
            anchor="picard-illuminalanemetrics",
            description="Quality control metrics on cluster density for each lane of an Illumina flowcell. For more information, see the [Picard Documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html#IlluminaLaneMetrics).",
            plot=lane_metrics_table(self.picard_lane_metrics),
        )

    # Return the number of detected samples to the parent module
    return len(self.picard_lane_metrics)
