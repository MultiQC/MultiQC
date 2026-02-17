"""MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging
from typing import Any, Dict

from multiqc.plots import linegraph

from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """
    Find Picard QualityByCycleMetrics reports and parse their data.

    Supports both standard format (CYCLE, MEAN_QUALITY) and extended format
    (CYCLE, MEAN_QUALITY, MEAN_ORIGINAL_QUALITY) for BQSR comparison.

    Returns:
        Set of sample names that were successfully parsed
    """

    # First try the extended format with MEAN_ORIGINAL_QUALITY
    headers_extended = ["CYCLE", "MEAN_QUALITY", "MEAN_ORIGINAL_QUALITY"]
    formats_extended = [int, float, float]
    all_data_extended = read_histogram(
        self,
        "picard/quality_by_cycle",
        headers_extended,
        formats_extended,
        picard_tool="MeanQualityByCycle",
        sentieon_algo="MeanQualityByCycle",
    )

    # Then try the standard format for any remaining files
    headers_standard = ["CYCLE", "MEAN_QUALITY"]
    formats_standard = [int, float]
    all_data_standard = read_histogram(
        self,
        "picard/quality_by_cycle",
        headers_standard,
        formats_standard,
        picard_tool="MeanQualityByCycle",
        sentieon_algo="MeanQualityByCycle",
    )

    # Combine all data - merge intelligently to avoid data loss
    all_data = {}
    has_original_quality = False

    # Start with standard format data
    all_data.update(all_data_standard)

    # Add extended format data, only overwriting if sample has extended data
    if all_data_extended:
        extended_samples_count = 0
        # Only overwrite if the extended data is actually extended (has MEAN_ORIGINAL_QUALITY)
        for s_name, s_data in all_data_extended.items():
            if s_data:  # Ensure sample has data
                sample_data: Dict[str, Any] = next(iter(s_data.values()), {})
                if "MEAN_ORIGINAL_QUALITY" in sample_data:
                    all_data[s_name] = s_data
                    has_original_quality = True
                    extended_samples_count += 1
                elif s_name not in all_data:
                    # If sample wasn't in standard format, still add it
                    all_data[s_name] = s_data

        if extended_samples_count > 0:
            log.debug(f"Found {extended_samples_count} samples with MEAN_ORIGINAL_QUALITY data for BQSR comparison")

    if not all_data:
        return set()

    # Filter out empty samples
    all_data = {s_name: s_data for s_name, s_data in all_data.items() if s_data}

    if not all_data:
        return set()

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
        "x_decimals": False,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f}",
        "ymin": 0,
    }

    # Prepare data for plotting
    lg_mean_quality = {}
    lg_original_quality = {}

    for s_name in all_data:
        if not all_data[s_name]:
            continue

        # All samples get mean quality data
        lg_mean_quality[s_name] = {cycle: data["MEAN_QUALITY"] for cycle, data in all_data[s_name].items()}

        # Only samples with MEAN_ORIGINAL_QUALITY get original quality data
        try:
            sample_data = next(iter(all_data[s_name].values()))
            if "MEAN_ORIGINAL_QUALITY" in sample_data:
                # Validate that all cycles have MEAN_ORIGINAL_QUALITY data
                original_data = {}
                for cycle, data in all_data[s_name].items():
                    if "MEAN_ORIGINAL_QUALITY" in data and data["MEAN_ORIGINAL_QUALITY"] is not None:
                        original_data[cycle] = data["MEAN_ORIGINAL_QUALITY"]

                if original_data:  # Only add if we have valid original quality data
                    lg_original_quality[s_name] = original_data
        except StopIteration:
            log.warning(f"No cycle data found for sample {s_name}")
            continue

    if has_original_quality and lg_original_quality:
        # Plot both mean quality and mean original quality
        plot_data = [lg_mean_quality, lg_original_quality]

        # Update plot config with data labels
        pconfig.update(
            {
                "data_labels": [
                    {"name": "Recalibrated Quality", "ylab": "Mean Base Quality"},
                    {"name": "Original Quality", "ylab": "Mean Base Quality"},
                ]
            }
        )

        description = "Plot shows the mean base quality by cycle, comparing recalibrated quality scores with original quality scores."
        helptext_addition = """
        
        When Base Quality Score Recalibration (BQSR) has been applied, this plot shows both the 
        recalibrated quality scores (Mean Quality) and the original quality scores (Mean Original Quality). 
        This comparison helps assess the effectiveness of BQSR or other quality adjustments.
        """
    else:
        # Plot only mean quality (legacy behavior)
        plot_data = [lg_mean_quality]
        description = "Plot shows the mean base quality by cycle."
        helptext_addition = ""

    self.add_section(
        name="Mean Base Quality by Cycle",
        anchor=f"{self.anchor}-quality-by-cycle",
        description=description,
        helptext=f"""
        This metric gives an overall snapshot of sequencing machine performance.
        For most types of sequencing data, the output is expected to show a slight
        reduction in overall base quality scores towards the end of each read.

        Spikes in quality within reads are not expected and may indicate that technical
        problems occurred during sequencing.{helptext_addition}
        """,
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Return the number of detected samples to the parent module
    return all_data.keys()
