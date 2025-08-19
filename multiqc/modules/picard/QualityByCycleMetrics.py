"""MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging

from multiqc.plots import linegraph

from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard QualityByCycleMetrics reports and parse their data"""

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

    # Combine all data - extended format takes precedence for samples that have it
    all_data = {}
    has_original_quality = False

    # Add standard format data first
    all_data.update(all_data_standard)

    # Add extended format data (may overwrite standard format for same samples)
    if all_data_extended:
        all_data.update(all_data_extended)
        has_original_quality = True

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
    if has_original_quality:
        # Plot both mean quality and mean original quality
        lg_mean_quality = {}
        lg_original_quality = {}

        for s_name in all_data:
            # All samples get mean quality data
            lg_mean_quality[s_name] = dict((cycle, data["MEAN_QUALITY"]) for cycle, data in all_data[s_name].items())

            # Only samples with MEAN_ORIGINAL_QUALITY get original quality data
            if "MEAN_ORIGINAL_QUALITY" in next(iter(all_data[s_name].values())):
                lg_original_quality[s_name] = dict(
                    (cycle, data["MEAN_ORIGINAL_QUALITY"]) for cycle, data in all_data[s_name].items()
                )

        # Combine both datasets with data labels
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
        lg_mean_quality = {}
        for s_name in all_data:
            lg_mean_quality[s_name] = dict((cycle, data["MEAN_QUALITY"]) for cycle, data in all_data[s_name].items())

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
