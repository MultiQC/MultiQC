""" MultiQC submodule to parse output from Picard CollectIlluminaBasecallingMetrics """

import logging
from collections import OrderedDict

from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard CollectIlluminaBasecallingMetrics reports and parse their data"""

    # Set up vars
    self.picard_basecalling_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/collectilluminabasecallingmetrics", filehandles=True):
        self.add_data_source(f, section="CollectIlluminaBasecallingMetrics")

        keys = None
        for line in f["f"]:
            if "IlluminaBasecallingMetrics" and "## METRICS CLASS" in line:
                keys = f["f"].readline().strip("\n").split("\t")
            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) == len(keys):
                    data = dict(zip(keys, vals))
                    if not data["MOLECULAR_BARCODE_SEQUENCE_1"]:
                        data.pop("MOLECULAR_BARCODE_SEQUENCE_1")
                        data.pop("MOLECULAR_BARCODE_NAME")
                        self.picard_basecalling_metrics[data["LANE"]] = data

    # Filter to strip out ignored sample names
    self.picard_basecalling_metrics = self.ignore_samples(self.picard_basecalling_metrics)

    if len(self.picard_basecalling_metrics) > 0:
        # Write parsed data to a file
        self.write_data_file(self.picard_basecalling_metrics, "multiqc_picard_IlluminaBasecallingMetrics")

        self.add_section(
            name="Basecalling Metrics",
            anchor="picard-illuminabasecallingmetrics",
            description="Quality control metrics for each lane of an Illumina flowcell.",
            helptext="""
            For full details, please see the [Picard Documentation](http://broadinstitute.github.io/picard/picard-metric-definitions.html#IlluminaBasecallingMetrics).

            * `PF_BASES` / `NPF_BASES` :  The total number of passing-filter / not-passing-filter bases assigned to the index.
            * `PF_READS` / `NPF_READS` :  The total number of passing-filter / not-passing-filter reads assigned to the index.
            * `PF_CLUSTERS` / `NPF_CLUSTERS` :  The total number of passing-filter / not-passing-filter clusters assigned to the index.

            `NPF` stands for _"not passing filter"_ and is calculated by subtracting the `PF_` metric from the `TOTAL_` Picard metrics.
            """,
            plot=lane_metrics_plot(self.picard_basecalling_metrics),
        )

    # Return the number of detected samples to the parent module
    return len(self.picard_basecalling_metrics)


def lane_metrics_table(data):
    headers = OrderedDict()
    headers["TOTAL_BASES"] = {"title": "Total Bases"}
    headers["PF_BASES"] = {"title": "Passing Filter Bases"}
    headers["TOTAL_READS"] = {"title": "Total Reads"}
    headers["PF_READS"] = {"title": "Passing Filter Reads"}
    headers["TOTAL_CLUSTERS"] = {"title": "Total Cluster"}
    headers["PF_CLUSTERS"] = {"title": "Passing Filter Clusters"}

    table_config = {
        "id": "picard-illumina-basecalling-metrics-table",
        "namespace": "Picard",
        "table_title": "Picard Illumina Base Calling Metrics",
    }
    tdata = {}
    for lane_number, lane in data.items():
        tdata[f"L{lane_number}"] = lane
    return table.plot(tdata, headers, table_config)


def lane_metrics_plot(data):
    plot_config = {
        "id": "plot-picard-illuminabasecallingmetrics",
        "title": "Picard: Illumina Basecalling Metrics",
        "ylab": "Lane",
        "data_labels": [
            {"name": "Bases", "ylab": "Number of Bases"},
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Clusters", "ylab": "Number of Clusters"},
        ],
    }

    plot_cats = [OrderedDict(), OrderedDict(), OrderedDict()]
    plot_cats[0]["PF_BASES"] = {"title": "Passing Filter Bases"}
    plot_cats[0]["NPF_BASES"] = {"title": "Non Passing Filter Bases"}
    plot_cats[1]["PF_READS"] = {"title": "Passing Filter Reads"}
    plot_cats[1]["NPF_READS"] = {"title": "Non Passing Filter Reads"}
    plot_cats[2]["PF_CLUSTERS"] = {"title": "Passing Filter Clusters"}
    plot_cats[2]["NPF_CLUSTERS"] = {"title": "Non Passing Filter Clusters"}
    tdata = {}
    for lane_number, lane in data.items():
        lane["NPF_BASES"] = int(lane["TOTAL_BASES"]) - int(lane["PF_BASES"])
        lane["NPF_READS"] = int(lane["TOTAL_READS"]) - int(lane["PF_READS"])
        lane["NPF_CLUSTERS"] = int(lane["TOTAL_CLUSTERS"]) - int(lane["PF_CLUSTERS"])
        tdata[f"L{lane_number}"] = lane
    return bargraph.plot([tdata, tdata, tdata], plot_cats, plot_config)
