""" MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging

from multiqc.plots import linegraph

from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard QualityByCycleMetrics reports and parse their data"""

    headers = ["CYCLE", "MEAN_QUALITY"]
    formats = [int, float]
    all_data = read_histogram(
        self,
        "picard/quality_by_cycle",
        headers,
        formats,
        picard_tool="MeanQualityByCycle",
        sentieon_algo="MeanQualityByCycle",
    )

    if not all_data:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Write parsed data to a file
    self.write_data_file(all_data, f"multiqc_{self.anchor}_quality_by_cycle")

    # Plot the data and add section
    pconfig = {
        "id": f"{self.anchor}_quality_by_cycle",
        "title": f"{self.name}: Mean Base Quality by Cycle",
        "ylab": "Mean Base Quality",
        "xlab": "Cycle Number",
        "xDecimals": False,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f}",
        "ymin": 0,
    }

    lg = {}
    for s_name in all_data:
        lg[s_name] = dict((cycle, data["MEAN_QUALITY"]) for cycle, data in all_data[s_name].items())

    self.add_section(
        name="Mean Base Quality by Cycle",
        anchor=f"{self.anchor}-quality-by-cycle",
        description="Plot shows the mean base quality by cycle.",
        helptext="""
        This metric gives an overall snapshot of sequencing machine performance.
        For most types of sequencing data, the output is expected to show a slight
        reduction in overall base quality scores towards the end of each read.

        Spikes in quality within reads are not expected and may indicate that technical
        problems occurred during sequencing.
        """,
        plot=linegraph.plot([lg], pconfig),
    )

    # Return the number of detected samples to the parent module
    return len(all_data)
