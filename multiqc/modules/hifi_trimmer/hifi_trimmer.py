import json
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.samtools.stats import parse_samtools_stats_lines
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    A command-line tool for filtering and trimming extraneous adapter hits from
    a HiFi read set using a BLAST search against a fasta file of adapter sequences.
    It is designed to be highly configurable, with per-adapter settings to determine
    actions if the adapter is found at the ends of a read or in the middle.
    A `samtools stats` file can be provided to give a more accurate count of the
    number of reads / bases involved for this sample.
    Linking is done by sample name.
    """

    def __init__(self):
        super().__init__(
            name="HiFi-Trimmer",
            anchor="hifi_trimmer",
            href="https://github.com/sanger-tol/hifi-trimmer",
            info="Filters and trims adapter sequences from HiFi reads using BLAST.",
            doi="",  # No DOI available
        )

        # Find and load any HiFi-Trimmer reports
        self.hifi_trimmer_data = {}

        # First, parse samtools stats files to get total sequences and bases
        samtools_data = {}
        for f in self.find_log_files("samtools/stats"):
            parsed_data, samtools_version, htslib_version = parse_samtools_stats_lines(f["f"])
            if "raw_total_sequences" in parsed_data or "total_length" in parsed_data:
                s_name = f["s_name"]
                samtools_data[s_name] = parsed_data
                log.debug(f"Found samtools stats for {s_name}")

        # Then parse HiFi-Trimmer JSON files
        for f in self.find_log_files("hifi_trimmer", filehandles=True):
            parsed = self.parse_hifi_trimmer_json(f)
            if parsed is not None:
                s_name = parsed["s_name"]
                if s_name in self.hifi_trimmer_data:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                data = parsed["data"]

                # Try to merge with samtools stats data if available
                # Try exact match first, then try without common suffixes
                samtools_key = None
                if s_name in samtools_data:
                    samtools_key = s_name
                else:
                    # Try stripping common suffixes to match samtools stats files
                    for suffix in ["_hifi_trimmer", ".hifi_trimmer", "_trimmer"]:
                        if s_name.endswith(suffix):
                            candidate = s_name[: -len(suffix)]
                            if candidate in samtools_data:
                                samtools_key = candidate
                                break

                if samtools_key:
                    # Use samtools stats data as the source of truth for total reads/bases
                    data["sample_total_reads"] = samtools_data[samtools_key].get("raw_total_sequences", 0)
                    data["sample_total_bases"] = samtools_data[samtools_key].get("total_length", 0)
                    log.debug(f"Merged samtools stats data for {s_name} (matched with {samtools_key})")

                self.hifi_trimmer_data[s_name] = data
                self.add_data_source(f, s_name=s_name)

        # Filter to strip out ignored sample names
        self.hifi_trimmer_data = self.ignore_samples(self.hifi_trimmer_data)

        if len(self.hifi_trimmer_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hifi_trimmer_data)} reports")

        # Register HiFi-Trimmer in the software versions table; passing None is a supported no-op
        # when the version cannot be determined from the available reports.
        self.add_software_version(None)

        # Add data to general stats table
        self.hifi_trimmer_general_stats_table()

        # Add sections with bar charts
        self.add_section(
            name="Read Statistics",
            anchor="hifi_trimmer_reads",
            description="Summary of read processing statistics from HiFi-Trimmer.",
            helptext="This plot shows the number of reads that were processed, trimmed, discarded by HiFi-Trimmer &mdash; and unprocessed, if total read counts are available from samtools stats.",
            plot=self.hifi_trimmer_reads_barplot(),
        )

        self.add_section(
            name="Base Statistics",
            anchor="hifi_trimmer_bases",
            description="Summary of base processing statistics from HiFi-Trimmer.",
            helptext="This plot shows the number of bases that were processed, kept, removed by HiFi-Trimmer &mdash; and unprocessed, if total base counts are available from samtools stats.",
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
            "sample_total_reads": {
                "title": "Total Reads",
                "description": "Total number of reads in the sample (from samtools stats)",
                "min": 0,
                "scale": "Greens",
                "format": "{:,.0f}",
                "shared_key": "read_count",
                "hidden": True,
            },
            "total_reads_processed": {
                "title": "Reads Processed",
                "description": "Total number of reads processed by HiFi-Trimmer",
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
            "sample_total_bases": {
                "title": "Total Bases",
                "description": "Total number of bases in the sample (from samtools stats)",
                "min": 0,
                "scale": "Purples",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "total_bases_processed": {
                "title": "Bases Processed",
                "description": "Total number of bases processed by HiFi-Trimmer",
                "min": 0,
                "scale": "Purples",
                "format": "{:,.0f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.hifi_trimmer_data, headers)

    def hifi_trimmer_reads_barplot(self):
        """Generate a horizontal bar plot showing read statistics"""

        # Prepare data with calculated unprocessed reads
        plot_data = {}
        for s_name, data in self.hifi_trimmer_data.items():
            plot_data[s_name] = {
                "total_reads_kept": data.get("total_reads_kept", 0),
                "total_reads_trimmed": data.get("total_reads_trimmed", 0),
                "total_reads_discarded": data.get("total_reads_discarded", 0),
            }
            # Add unprocessed reads if we have sample totals from samtools stats
            if "sample_total_reads" in data:
                sample_total = data["sample_total_reads"]
                processed_total = data.get("total_reads_processed", 0)
                if sample_total > processed_total:
                    plot_data[s_name]["unprocessed_reads"] = sample_total - processed_total

        cats = {
            "unprocessed_reads": {"name": "Unprocessed Reads", "color": "#95a5a6"},
            "total_reads_kept": {"name": "Reads Kept", "color": "#2ecc71"},
            "total_reads_trimmed": {"name": "Reads Trimmed", "color": "#f39c12"},
            "total_reads_discarded": {"name": "Reads Discarded", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_reads_plot",
            "title": "HiFi-Trimmer: Read Statistics",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_zero_cats": True,
        }

        return bargraph.plot(plot_data, cats, config)

    def hifi_trimmer_bases_barplot(self):
        """Generate a horizontal bar plot showing base statistics"""

        # Prepare data with calculated unprocessed bases
        plot_data = {}
        for s_name, data in self.hifi_trimmer_data.items():
            plot_data[s_name] = {
                "total_bases_kept": data.get("total_bases_kept", 0),
                "total_bases_removed": data.get("total_bases_removed", 0),
            }
            # Add unprocessed bases if we have sample totals from samtools stats
            if "sample_total_bases" in data:
                sample_total = data["sample_total_bases"]
                processed_total = data.get("total_bases_processed", 0)
                if sample_total > processed_total:
                    plot_data[s_name]["unprocessed_bases"] = sample_total - processed_total

        cats = {
            "unprocessed_bases": {"name": "Unprocessed Bases", "color": "#95a5a6"},
            "total_bases_kept": {"name": "Bases Kept", "color": "#3498db"},
            "total_bases_removed": {"name": "Bases Removed", "color": "#e74c3c"},
        }

        config = {
            "id": "hifi_trimmer_bases_plot",
            "title": "HiFi-Trimmer: Base Statistics",
            "ylab": "Number of Bases",
            "cpswitch_counts_label": "Number of Bases",
            "hide_zero_cats": True,
        }

        return bargraph.plot(plot_data, cats, config)
