"""MultiQC submodule to parse output from Picard QualityScoreDistribution"""

import logging
from collections import OrderedDict

from multiqc.plots import linegraph

from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard QualityScoreDistribution reports and parse their data"""

    headers = ["QUALITY", "COUNT_OF_Q"]
    formats = [int, int]
    all_data = read_histogram(
        self,
        "picard/quality_score_distribution",
        headers,
        formats,
        picard_tool="QualityScoreDistribution",
        sentieon_algo="QualDistribution",
    )

    if not all_data:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Write parsed data to a file
    self.write_data_file(all_data, f"multiqc_{self.anchor}_quality_score_distribution")

    # Plot the data and add section
    pconfig = {
        "id": f"{self.anchor}_quality_score_distribution",
        "title": f"{self.name}: Base Quality Distribution",
        "ylab": "Number of Bases",
        "xlab": "Base Quality Score",
        "x_decimals": False,
        "tt_label": "<b>base quality{point.x}</b>: {point.y}",
        "ymin": 0,
    }

    lg = {}
    for s_name in all_data:
        lg[s_name] = OrderedDict((qual, data["COUNT_OF_Q"]) for qual, data in all_data[s_name].items())

    self.add_section(
        name="Base Quality Distribution",
        anchor=f"{self.anchor}-quality-score-distribution",
        description="Plot shows the count of each base quality score.",
        plot=linegraph.plot([lg], pconfig),
    )

    # Return the number of detected samples to the parent module
    return all_data.keys()
