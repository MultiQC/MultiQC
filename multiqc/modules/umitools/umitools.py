"""MultiQC module to parse output from UMI-tools"""

import logging
import re
from typing import Dict, Optional

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    umitools module class, parses dedup logs
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="UMI-tools",
            anchor="umitools",
            href="https://github.com/CGATOxford/UMI-tools",
            info="contains tools for dealing with Unique Molecular Identifiers (UMIs)/(RMTs) and scRNA-Seq barcodes.",
            doi="10.1101/gr.209601.116",
        )

        dedup_data_by_sample = dict()
        for f in self.find_log_files("umitools/dedup"):
            # Parse the log file for sample name and statistics
            data = self.parse_dedup_logs(f)
            if data:
                f["s_name"] = self._parse_s_name(f) or f["s_name"]
                if f["s_name"] in dedup_data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                dedup_data_by_sample[f["s_name"]] = data
                self.add_data_source(f)  # Add the log file information to the multiqc_sources.txt
                if "version" in data:  # Add version info
                    self.add_software_version(data.pop("version"), f["s_name"])

        extract_data_by_sample = dict()
        for f in self.find_log_files("umitools/extract"):
            # Parse the log file for sample name and statistics
            data = self.parse_extract_logs(f)
            if data:
                f["s_name"] = self._parse_s_name(f) or f["s_name"]
                if f["s_name"] in extract_data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                extract_data_by_sample[f["s_name"]] = data
                self.add_data_source(f)  # Add the log file information to the multiqc_sources.txt
                if "version" in data:  # Add version info
                    self.add_software_version(data.pop("version"), f["s_name"])

        # Check the log files against the user supplied list of samples to ignore
        dedup_data_by_sample = self.ignore_samples(dedup_data_by_sample)
        extract_data_by_sample = self.ignore_samples(extract_data_by_sample)
        if len(dedup_data_by_sample) == 0 and len(extract_data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(dedup_data_by_sample)} deduplication reports")
        log.info(f"Found {len(extract_data_by_sample)} extract reports")

        # Write parsed reports data to a file
        self.write_data_file(dedup_data_by_sample, "multiqc_umitools_dedup")
        self.write_data_file(extract_data_by_sample, "multiqc_umitools_extract")

        if dedup_data_by_sample:
            self.dedup_general_stats_table(dedup_data_by_sample)
            self.deduplication_plot(dedup_data_by_sample)
            self.umi_stats_violin(dedup_data_by_sample)

        if extract_data_by_sample:
            self.extract_general_stats_table(extract_data_by_sample)
            self.umitools_extract_barplot(extract_data_by_sample)

            # Only plot regex match rate if we have regex extraction logs:
            if any(
                any(key in v for key in ["read1_match", "read2_match", "read1_mismatch", "read2_mismatch"])
                for v in extract_data_by_sample.values()
            ):
                self.umitools_extract_barplot_regex(extract_data_by_sample)

    def _parse_s_name(self, f) -> Optional[str]:
        # Get the s_name from the input file if possible
        # stdin : <_io.TextIOWrapper name='M18-39155_T1.Aligned.sortedByCoord.out.bam' mode='r' encoding='UTF-8'>
        s_name_re = r"stdin\s+:\s+<_io\.TextIOWrapper name='([^\']+)'"
        s_name_match = re.search(s_name_re, f["f"])
        if s_name_match:
            return self.clean_s_name(s_name_match.group(1))
        return None

    @staticmethod
    def parse_dedup_logs(f) -> Dict:
        regexes = [
            (int, "total_umis", r"INFO total_umis (\d+)"),
            (int, "unique_umis", r"INFO #umis (\d+)"),
            (int, "dedup_input_reads", r"INFO Reads: Input Reads: (\d+)"),
            (int, "dedup_output_reads", r"INFO Number of reads out: (\d+)"),
            (int, "positions_deduplicated", r"INFO Total number of positions deduplicated: (\d+)"),
            (float, "mean_umi_per_pos", r"INFO Mean number of unique UMIs per position: ([\d\.]+)"),
            (float, "max_umi_per_pos", r"INFO Max. number of unique UMIs per position: (\d+)"),
            (str, "version", r"# UMI-tools version: ([\d\.]+)"),
        ]

        data = {}
        # Search for values using regular expressions
        for type_, key, regex in regexes:
            re_matches = re.search(regex, f["f"])
            if re_matches:
                if key == "version":
                    data[key] = re_matches.group(1)
                else:
                    data[key] = type_(re_matches.group(1))

            # Calculate a few simple supplementary stats
            try:
                data["dedup_percent_passing"] = round(
                    ((data["dedup_output_reads"] / data["dedup_input_reads"]) * 100.0), 2
                )
            except (KeyError, ZeroDivisionError):
                pass
            try:
                data["dedup_removed_reads"] = data["dedup_input_reads"] - data["dedup_output_reads"]
            except KeyError:
                pass

        return data

    @staticmethod
    def parse_extract_logs(f) -> Dict:
        regexes = [
            (int, "extract_input_reads", r"INFO Input Reads: (\d+)"),
            (int, "read1_mismatch", r"INFO regex does not match read1: (\d+)"),
            (int, "read1_match", r"INFO regex matches read1: (\d+)"),
            (int, "read2_mismatch", r"INFO regex does not match read2: (\d+)"),
            (int, "read2_match", r"INFO regex matches read2: (\d+)"),
            (int, "extract_output_reads", r"INFO Reads output: (\d+)"),
            (str, "version", r"# UMI-tools version: ([\d\.]+)"),
        ]

        data = {}
        # Search for values using regular expressions
        for type_, key, regex in regexes:
            re_matches = re.search(regex, f["f"])
            if re_matches:
                if key == "version":
                    data[key] = re_matches.group(1)
                else:
                    data[key] = type_(re_matches.group(1))

            # Calculate a few simple supplementary stats
            try:
                data["extract_percent_passing"] = round(
                    ((data["extract_output_reads"] / data["extract_input_reads"]) * 100.0), 2
                )
            except (KeyError, ZeroDivisionError):
                pass
            try:
                data["extract_removed_reads"] = data["extract_input_reads"] - data["extract_output_reads"]
            except KeyError:
                pass

        return data

    def dedup_general_stats_table(self, data_by_sample):
        """
        Take the parsed stats from the umitools report and add it to the
        basic stats table at the top of the report
        """
        headers = {
            "dedup_output_reads": {
                "title": f"{config.read_count_prefix} Unique Reads",
                "description": f"Reads remaining after deduplication ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "PuRd",
            },
            "dedup_percent_passing": {
                "title": "% Pass Dedup",
                "description": "% processed reads that passed deduplication",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def extract_general_stats_table(self, data_by_sample):
        headers = {
            "extract_output_reads": {
                "title": f"{config.read_count_prefix} Extracted Reads",
                "description": f"Reads remaining after extraction ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "PuRd",
            },
            "extract_percent_passing": {
                "title": "% UMI Extract",
                "description": "% reads from which a UMI extraction succeeded.",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def deduplication_plot(self, data_by_sample):
        """Generate a plot the read deduplication rates for the main report"""

        # Specify the order of the different possible categories
        keys = {
            "dedup_output_reads": {"color": "#7fc97f", "name": "Reads remaining"},
            "dedup_removed_reads": {"color": "#fdc086", "name": "Reads removed"},
        }

        return self.add_section(
            name="Deduplicated Reads",
            anchor="umitools-dedup-plot",
            plot=bargraph.plot(
                data_by_sample,
                keys,
                {
                    "id": "umitools_deduplication_barplot",
                    "title": "UMI-tools: Deduplication Counts",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of Reads",
                },
            ),
        )

    def umitools_extract_barplot(self, data_by_sample):
        keys = {
            "extract_output_reads": {"color": "#7fc9c7", "name": "Extraction succeeded"},
            "extract_removed_reads": {"color": "#fd9286", "name": "Extraction failed"},
        }

        # Add a section with a barplot plot of UMI stats to the report
        self.add_section(
            name="UMI extraction rate",
            anchor="umitools_extract_overall",
            description="Success rate of `umi_tools extract`",
            plot=bargraph.plot(
                data_by_sample,
                keys,
                {
                    "id": "umitools_extract_barplot_success",
                    "title": "UMI-tools: Extraction rate",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of Reads",
                },
            ),
        )

    def umitools_extract_barplot_regex(self, data_by_sample):
        keys = [
            {
                "read1_match": {"color": "#7fc9c7", "name": "Match on read 1"},
                "read1_mismatch": {"color": "#fd9286", "name": "Mismatch on read 1"},
            },
            {
                "read2_match": {"color": "#7fc9c7", "name": "Match on read 2"},
                "read2_mismatch": {"color": "#fd9286", "name": "Mismatch on read 2"},
            },
        ]

        subset_read1 = {
            sample: {key: metrics[key] for key in ("read1_match", "read1_mismatch") if key in metrics}
            for sample, metrics in data_by_sample.items()
        }

        subset_read2 = {
            sample: {key: metrics[key] for key in ("read2_match", "read2_mismatch") if key in metrics}
            for sample, metrics in data_by_sample.items()
        }

        # Add a section with a barplot plot of the regex match rate to the report
        self.add_section(
            name="UMI extraction regex match rate",
            anchor="umitools_extract_regex",
            description="Regex match rate of `umi_tools extract`",
            plot=bargraph.plot(
                [subset_read1, subset_read2],
                keys,
                {
                    "id": "umitools_extract_regex_barplot",
                    "title": "UMI-tools: Extraction regex match rate",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of Reads",
                    "data_labels": ["Read 1", "Read 2"],
                },
            ),
        )

    def umi_stats_violin(self, data_by_sample):
        """Generate a violin of UMI stats for the main report"""

        headers = {
            "positions_deduplicated": {
                "title": "Positions Dedup",
                "description": "Total number of positions deduplicated",
                "min": 0,
                "format": "{:,d}",
                "scale": "Greens",
            },
            "total_umis": {
                "title": "Total UMIs",
                "description": "Total UMIs found in sample",
                "min": 0,
                "format": "{:,d}",
                "scale": "Blues",
            },
            "unique_umis": {
                "title": "Unique UMIs",
                "description": "Unique UMIs found in sample",
                "min": 0,
                "format": "{:,d}",
                "scale": "Purples",
            },
            "mean_umi_per_pos": {
                "title": "Mean #UMI",
                "description": "Mean number of unique UMIs per position",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Reds",
            },
            "max_umi_per_pos": {
                "title": "Max #UMI",
                "description": "Max number of unique UMIs per position",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
        }

        # add a section with a violin plot of UMI stats to the report
        self.add_section(
            name="UMI Stats",
            anchor="umitools-umi-stats",
            description="Statistics from running `umi_tools dedup` or `umi_tools extract`",
            helptext="""
            - **Positions Dedup**: Total number of positions deduplicated
            - **Total UMIs**: Total UMIs found in sample
            - **Unique UMIs**: Unique UMIs found in sample
            - **Mean #UMI**: Mean number of unique UMIs per position
            - **Max #UMI**: Max number of unique UMIs per position
            """,
            plot=violin.plot(
                data_by_sample,
                headers,
                {
                    "id": "umitools_stats_violin",
                    "title": "UMI-tools: UMI stats",
                },
            ),
        )
