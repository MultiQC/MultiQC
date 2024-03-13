""" MultiQC submodule to parse output from Picard IlluminaBasecallingMetrics """

import logging

from multiqc.modules.picard import util
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard IlluminaBasecallingMetrics reports and parse their data"""

    data_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/collectilluminabasecallingmetrics", filehandles=True):
        keys = None

        for line in f["f"]:
            if util.is_line_right_before_table(line, "IlluminaBasecallingMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys):
                    keys = None
                    continue

                data = dict(zip(keys, vals))
                # We only care about the last line
                if data["MOLECULAR_BARCODE_SEQUENCE_1"].strip() == "":
                    data.pop("MOLECULAR_BARCODE_SEQUENCE_1")
                    data.pop("MOLECULAR_BARCODE_NAME")
                    s_name = data["LANE"]
                    if s_name in data_by_sample:
                        log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                    module.add_data_source(f, s_name=s_name, section="IlluminaBasecallingMetrics")
                    data_by_sample[s_name] = data

    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_IlluminaBasecallingMetrics")

    module.add_section(
        name="Basecalling Metrics",
        anchor=f"{module.anchor}-illuminabasecallingmetrics",
        description="Quality control metrics for each lane of an Illumina flowcell.",
        helptext="""
        For full details, please see the [Picard Documentation](http://broadinstitute.github.io/picard/picard-metric-definitions.html#IlluminaBasecallingMetrics).

        * `PF_BASES` / `NPF_BASES` :  The total number of passing-filter / not-passing-filter bases assigned to the index.
        * `PF_READS` / `NPF_READS` :  The total number of passing-filter / not-passing-filter reads assigned to the index.
        * `PF_CLUSTERS` / `NPF_CLUSTERS` :  The total number of passing-filter / not-passing-filter clusters assigned to the index.

        `NPF` stands for _"not passing filter"_ and is calculated by subtracting the `PF_` metric from the `TOTAL_` Picard metrics.
        """,
        plot=lane_metrics_plot(module, data_by_sample),
    )

    # Return the number of detected samples to the parent module
    return len(data_by_sample)


def lane_metrics_table(self, data):
    headers = {
        "TOTAL_BASES": {"title": "Total Bases"},
        "PF_BASES": {"title": "Passing Filter Bases"},
        "TOTAL_READS": {"title": "Total Reads"},
        "PF_READS": {"title": "Passing Filter Reads"},
        "TOTAL_CLUSTERS": {"title": "Total Cluster"},
        "PF_CLUSTERS": {"title": "Passing Filter Clusters"},
    }

    table_config = {
        "id": f"{self.anchor}-illumina-basecalling-metrics-table",
        "namespace": self.name,
        "table_title": f"{self.name} Illumina Base Calling Metrics",
    }
    tdata = {}
    for lane_number, lane in data.items():
        tdata[f"L{lane_number}"] = lane
    return table.plot(tdata, headers, table_config)


def lane_metrics_plot(self, data):
    plot_config = {
        "id": f"plot-{self.anchor}-illuminabasecallingmetrics",
        "title": f"{self.name}: Illumina Basecalling Metrics",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Bases", "ylab": "Number of Bases"},
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Clusters", "ylab": "Number of Clusters"},
        ],
    }

    plot_cats = [
        {
            "PF_BASES": {"title": "Passing Filter Bases"},
            "NPF_BASES": {"title": "Non Passing Filter Bases"},
        },
        {
            "PF_READS": {"title": "Passing Filter Reads"},
            "NPF_READS": {"title": "Non Passing Filter Reads"},
        },
        {
            "PF_CLUSTERS": {"title": "Passing Filter Clusters"},
            "NPF_CLUSTERS": {"title": "Non Passing Filter Clusters"},
        },
    ]
    tdata = {}
    for lane_number, lane in data.items():
        lane["NPF_BASES"] = int(lane["TOTAL_BASES"]) - int(lane["PF_BASES"])
        lane["NPF_READS"] = int(lane["TOTAL_READS"]) - int(lane["PF_READS"])
        lane["NPF_CLUSTERS"] = int(lane["TOTAL_CLUSTERS"]) - int(lane["PF_CLUSTERS"])
        tdata[f"L{lane_number}"] = lane
    return bargraph.plot([tdata, tdata, tdata], plot_cats, plot_config)
