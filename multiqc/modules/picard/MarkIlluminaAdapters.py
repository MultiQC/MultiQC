""" MultiQC submodule to parse output from Picard MarkIlluminaAdapters """

import logging

from multiqc.plots import linegraph

from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard MarkIlluminaAdapters reports and parse their data"""

    headers = ["clipped_bases", "read_count"]
    formats = [int, int]
    data_by_sample = read_histogram(
        self,
        program_key="picard/markilluminaadapters",
        headers=headers,
        formats=formats,
        picard_tool="MarkIlluminaAdapters",
    )

    # Filter to strip out ignored sample names
    data_by_sample = self.ignore_samples(data_by_sample)
    if not data_by_sample:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Write parsed data to a file
    self.write_data_file(data_by_sample, "multiqc_picard_mark_illumina_adapters")

    # Plot the data and add section
    pconfig = {
        "id": "picard_mark_illumina_adapters",
        "title": "Picard: Mark Illumina Adapters",
        "ylab": "Clipped Bases",
        "xlab": "Cycle Number",
        "xDecimals": False,
        "tt_label": "<b>Cycle {point.x}</b>: {point.y:.2f}",
        "ymin": 0,
    }

    lg = {}
    for s_name in data_by_sample:
        lg[s_name] = {clipped_bases: data["read_count"] for clipped_bases, data in data_by_sample[s_name].items()}

    self.add_section(
        name="Mark Illumina Adapters",
        description="""
            Number of Clipped Bases by Read.
            See the [Picard Docuementation](https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters) for details.
        """,
        anchor="picard-mark-illumina-adapters",
        plot=linegraph.plot(lg, pconfig),
    )

    # Return the number of detected samples to the parent module
    return len(data_by_sample)
