"""MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging
from multiqc.plots import linegraph
from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """Find Picard QualityByCycleMetrics reports and parse their data"""

    # Attempt to read both MEAN_QUALITY and MEAN_ORIGINAL_QUALITY
    headers = ["CYCLE", "MEAN_QUALITY", "MEAN_ORIGINAL_QUALITY"]
    formats = [int, float, float]
    all_data = read_histogram(
        self,
        "picard/quality_by_cycle",
        headers,
        formats,
        picard_tool="MeanQualityByCycle",
        sentieon_algo="MeanQualityByCycle",
    )

    if not all_data:
        # Fallback to only MEAN_QUALITY if original quality isn't present
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
            return set()

        self.add_software_version(None)
        self.write_data_file(all_data, f"multiqc_{self.anchor}_quality_by_cycle")

        # Fallback: Only one line (mean quality)
        pconfig = {
            "id": f"{self.anchor}_quality_by_cycle",
            "title": f"{self.name}: Mean Base Quality by Cycle",
            "ylab": "Mean Base Quality",
            "xlab": "Cycle Number",
            "x_decimals": False,
            "tt_label": "<b>Cycle {point.x}</b>: {point.y:.2f}",
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

        return all_data.keys()

    # Proceed with full plot if both columns exist
    self.add_software_version(None)
    self.write_data_file(all_data, f"multiqc_{self.anchor}_quality_by_cycle")

    pconfig = {
        "id": f"{self.anchor}_quality_by_cycle",
        "title": f"{self.name}: Mean Base Quality by Cycle",
        "ylab": "Mean Base Quality",
        "xlab": "Cycle Number",
        "x_decimals": False,
        "tt_label": "<b>Cycle {point.x}</b>: {point.y:.2f}",
        "ymin": 0,
        "data_labels": [
            {"name": "Mean quality", "ylab": "Mean Base Quality"},
            {"name": "Original mean quality", "ylab": "Original Mean Base Quality"}
        ]
    }

    # Two line plots: one for each quality metric
    lg_mean = {}
    lg_orig = {}
    for s_name in all_data:
        lg_mean[s_name] = dict((cycle, data["MEAN_QUALITY"]) for cycle, data in all_data[s_name].items())
        lg_orig[s_name] = dict((cycle, data["MEAN_ORIGINAL_QUALITY"]) for cycle, data in all_data[s_name].items())

    self.add_section(
        name="Mean Base Quality by Cycle",
        anchor=f"{self.anchor}-quality-by-cycle",
        description="Plot shows both the mean and original mean base quality by sequencing cycle.",
        helptext="""
        This metric provides a visual overview of sequencing quality across cycles.
        The 'original' quality values (pre-recalibration or correction) help compare
        base calling confidence before and after adjustments.
        """,
        plot=linegraph.plot([lg_mean, lg_orig], pconfig),
    )

    return all_data.keys()
