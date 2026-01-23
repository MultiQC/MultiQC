import json
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The HiFi-Trimmer module parses JSON output files containing trimming statistics.
    It creates two horizontal bar charts: one showing read statistics (processed, discarded, trimmed)
    and another showing base statistics (processed, removed).
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HiFi-Trimmer",
            anchor="hifi_trimmer",
            href="https://github.com/treehouse-bioinformatics/hifi-trimmer",
            info="Trims adapters and low-quality bases from PacBio HiFi reads.",
        )

        # Find and load any HiFi-Trimmer reports
        self.hifi_trimmer_data = {}

        for f in self.find_log_files("hifi_trimmer", filehandles=True):
            parsed = self.parse_hifi_trimmer_json(f)
            if parsed is not None:
                s_name = parsed["s_name"]
                if s_name in self.hifi_trimmer_data:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                self.hifi_trimmer_data[s_name] = parsed["data"]
                self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names
        self.hifi_trimmer_data = self.ignore_samples(self.hifi_trimmer_data)

        if len(self.hifi_trimmer_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hifi_trimmer_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        self.add_software_version(None)

        # Add data to general stats table
        self.hifi_trimmer_general_stats_table()

        # Add sections with bar charts
        self.add_section(
            name="Read Statistics",
            anchor="hifi_trimmer_reads",
            description="Summary of read processing statistics from HiFi-Trimmer.",
            plot=self.hifi_trimmer_reads_barplot(),
        )

        self.add_section(
            name="Base Statistics",
            anchor="hifi_trimmer_bases",
            description="Summary of base processing statistics from HiFi-Trimmer.",
            plot=self.hifi_trimmer_bases_barplot(),
        )

        # Write parsed report data to a file
        self.write_data_file(self.hifi_trimmer_data, "multiqc_hifi_trimmer")

    def parse_hifi_trimmer_json(self, f):
        """Parse the JSON output from HiFi-Trimmer and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])

            # Check for required keys
            if "summary" not in parsed_json:
                log.debug(f"HiFi-Trimmer JSON missing 'summary' key - skipping sample: '{f['fn']}'")
                return None

            summary = parsed_json["summary"]
            required_keys = [
                "total_reads_processed",
                "total_bases_processed",
                "total_reads_discarded",
                "total_reads_trimmed",
                "total_bases_removed",
            ]

            if not all(key in summary for key in required_keys):
                log.debug(f"HiFi-Trimmer JSON missing required keys - skipping sample: '{f['fn']}'")
                return None

        except (json.JSONDecodeError, KeyError) as e:
            log.debug(f"Could not parse HiFi-Trimmer JSON: '{f['fn']}'")
            log.debug(e)
            return None

        # Get sample name from filename
        s_name = self.clean_s_name(f["fn"], f)

        # Extract metrics
        data = {
            "total_reads_processed": summary["total_reads_processed"],
            "total_bases_processed": summary["total_bases_processed"],
            "total_reads_discarded": summary["total_reads_discarded"],
            "total_reads_trimmed": summary["total_reads_trimmed"],
            "total_bases_removed": summary["total_bases_removed"],
        }

        # Calculate derived metrics
        data["total_reads_kept"] = data["total_reads_processed"] - data["total_reads_discarded"]
        data["total_bases_kept"] = data["total_bases_processed"] - data["total_bases_removed"]

        # Calculate percentages
        if data["total_reads_processed"] > 0:
            data["pct_reads_discarded"] = (data["total_reads_discarded"] / data["total_reads_processed"]) * 100
            data["pct_reads_trimmed"] = (data["total_reads_trimmed"] / data["total_reads_processed"]) * 100
            data["pct_reads_kept"] = (data["total_reads_kept"] / data["total_reads_processed"]) * 100
        else:
            data["pct_reads_discarded"] = 0
            data["pct_reads_trimmed"] = 0
            data["pct_reads_kept"] = 0

        if data["total_bases_processed"] > 0:
            data["pct_bases_removed"] = (data["total_bases_removed"] / data["total_bases_processed"]) * 100
            data["pct_bases_kept"] = (data["total_bases_kept"] / data["total_bases_processed"]) * 100
        else:
            data["pct_bases_removed"] = 0
            data["pct_bases_kept"] = 0

        return {"s_name": s_name, "data": data}

    def hifi_trimmer_general_stats_table(self):
        """Add key statistics to the general stats table"""
        headers = {
            "total_reads_processed": {
                "title": "Reads Processed",
                "description": "Total number of reads processed",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
                "shared_key": "read_count",
            },
            "pct_reads_kept": {
                "title": "% Reads Kept",
                "description": "Percentage of reads kept after filtering",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
            "total_bases_processed": {
                "title": "Bases Processed",
                "description": "Total number of bases processed",
                "min": 0,
                "scale": "Purples",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.hifi_trimmer_data, headers)

    def hifi_trimmer_reads_barplot(self):
        """Generate a horizontal bar plot showing read statistics"""

        cats = {
            "total_reads_kept": {"name": "Reads Kept", "color": "#2ecc71"},
            "total_reads_trimmed": {"name": "Reads Trimmed", "color": "#f39c12"},
            "total_reads_discarded": {"name": "Reads Discarded", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_reads_plot",
            "title": "HiFi-Trimmer: Read Statistics",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_zero_cats": False,
        }

        return bargraph.plot(self.hifi_trimmer_data, cats, config)

    def hifi_trimmer_bases_barplot(self):
        """Generate a horizontal bar plot showing base statistics"""

        cats = {
            "total_bases_kept": {"name": "Bases Kept", "color": "#3498db"},
            "total_bases_removed": {"name": "Bases Removed", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_bases_plot",
            "title": "HiFi-Trimmer: Base Statistics",
            "ylab": "Number of Bases",
            "cpswitch_counts_label": "Number of Bases",
            "hide_zero_cats": False,
        }

        return bargraph.plot(self.hifi_trimmer_data, cats, config)
