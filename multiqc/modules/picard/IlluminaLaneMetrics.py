""" MultiQC submodule to parse output from Picard CollectIlluminaLaneMetrics """

import logging
from collections import defaultdict

from multiqc.modules.picard import util
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


def lane_metrics_table(module, data):
    headers = {
        "CLUSTER_DENSITY": {
            "title": "Cluster Density",
            "description": "The number of clusters per unit area on this lane (cluster / " "mm^2`)",
            "scale": "Greens",
        },
        "TYPE_NAME": {
            "title": "Read Number",
            "description": "Defines an Illumina template read number (first or second)",
            "modify": lambda x: x.lower(),
        },
        "PREPHASING_APPLIED": {
            "title": "Prephasing Applied",
            "description": "Median pre-phasing value across all tiles in a lane, applied "
            "to the first and second template reads",
            "scale": "BuPu",
            "max": 1,
        },
        "PHASING_APPLIED": {
            "title": "Phasing Applied",
            "description": "Median phasing value across all tiles in a lane, applied to "
            "the first and second template reads",
            "scale": "BuPu",
            "max": 1,
        },
    }

    table_config = {
        "id": f"{module.anchor}-illumina-lane-metrics-table",
        "namespace": module.name,
        "table_title": f"{module.name} Illumina Lane Metrics",
    }
    tdata = {}
    for run_name, run in data.items():
        for lane_number, lane in run.items():
            tdata[f"{run_name} - L{lane_number}"] = lane
    return table.plot(tdata, headers, table_config)


def parse_reports(module):
    """Find Picard IlluminaLaneMetrics reports and parse their data"""

    # There can be two types of these files for the same sample: one with IlluminaLaneMetrics,
    # and one with IlluminaPhasingMetrics. We want to collect both, thus using default dicts
    # and calling .update() on them when any metrics are found.
    data_by_lane_by_run = defaultdict(lambda: defaultdict(dict))

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/collectilluminalanemetrics", filehandles=True):
        # Sample name from input file name by default
        run_name = f["s_name"]
        keys = None

        for line in f["f"]:
            maybe_run_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectIlluminaLaneMetrics",
                picard_opt="OUTPUT_PREFIX",
            )
            if maybe_run_name:
                run_name = maybe_run_name
                keys = None

            if run_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class=["IlluminaLaneMetrics", "IlluminaPhasingMetrics"]):
                keys = f["f"].readline().strip("\n").split("\t")
                if run_name in data_by_lane_by_run:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {run_name}")
                module.add_data_source(f, s_name=run_name, section="IlluminaLaneMetrics")

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys):
                    keys = None
                    continue

                d = dict(zip(keys, vals))
                lane = d["LANE"]
                data_by_lane_by_run[run_name][lane].update(d)

    data_by_lane_by_run = module.ignore_samples(data_by_lane_by_run)
    if len(data_by_lane_by_run) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_lane_by_run, f"multiqc_{module.anchor}_IlluminaLaneMetrics")

    module.add_section(
        name="Lane Metrics",
        anchor="picard-illuminalanemetrics",
        description="Quality control metrics on cluster density for each lane of an "
        "Illumina flowcell. For more information, see the [Picard "
        "Documentation]("
        "https://broadinstitute.github.io/picard/picard-metric"
        "-definitions.html#IlluminaLaneMetrics).",
        plot=lane_metrics_table(module, data_by_lane_by_run),
    )

    # Return the number of detected samples to the parent module
    return len(data_by_lane_by_run)
