"""MultiQC module to parse output from nanoq"""

import re
import logging
from typing import Dict, List, Any

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import mqc_colour

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Nanoq module"""

    _KEYS_MAPPING = {
        "Number of reads": "Number of reads",
        "Number of bases": "Number of bases",
        "N50 read length": "N50 read length",
        "Longest read": "Longest read",
        "Shortest read": "Shortest read",
        "Mean read length": "Mean read length",
        "Median read length": "Median read length",
        "Mean read quality": "Mean read quality",
        "Median read quality": "Median read quality",
    }

    _stat_types = ("summary", "quality")

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="nanoq",
            anchor="nanoq",
            href="https://github.com/nerdna/nanoq/",
            info="reports read quality and length from nanopore sequencing data in fast{a,q}.{gz,bz2,xz} format.",
            doi="10.21105/joss.02991",
        )

        # Find and load any nanoq reports
        self.nanoq_summary: Dict[str, Dict[str, float]] = {}
        self.nanoq_data: Dict[str, Dict[str, float]] = {}
        for f in self.find_log_files("nanoq", filehandles=True):
            self.parse_nanoq_log(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        # Filter to strip out ignored sample names
        self.nanoq_data = self.ignore_samples(self.nanoq_data)

        if len(self.nanoq_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.nanoq_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.nanoq_data, "multiqc_nanoq")

        # Add Nanoq summary to the general stats table
        self.nanoq_general_stats()

        # Quality distribution Plot
        self.reads_by_quality_plot()

        # Read length distribution Plot
        self.reads_by_length_plot()

    def parse_nanoq_log(self, f) -> None:
        """Parse output from nanoq"""
        # Initialize dictionaries to store the parsed data
        nanoq_summary: Dict[str, List[Any]] = {}
        read_length_thresholds: Dict[str, List[Any]] = {}
        read_quality_thresholds: Dict[str, List[Any]] = {}

        # Helper function to parse summary part
        def parse_summary(lines: List[str]) -> Dict[str, List[Any]]:
            summary: Dict[str, List[Any]] = {"Metric": [], "Value": []}
            for line in lines:
                if ":" in line:
                    metric, value = line.split(":", 1)
                    summary["Metric"].append(metric.strip())
                    summary["Value"].append(float(value.strip()))
            return summary

        # Helper function to parse thresholds part
        def parse_thresholds(lines: List[str], threshold_type: str) -> Dict[str, List[Any]]:
            thresholds: Dict[str, List[Any]] = {"Threshold": [], "Number of Reads": [], "Percentage": []}
            for line in lines:
                match = re.match(r">\s*(\d+)\s+(\d+)\s+(\d+\.\d+)%", line)
                if match:
                    threshold, num_reads, percentage = match.groups()
                    thresholds["Threshold"].append(threshold + (threshold_type if threshold_type == "bp" else ""))
                    thresholds["Number of Reads"].append(int(num_reads))
                    thresholds["Percentage"].append(float(percentage))
            return thresholds

        # Parse the file content
        segment = None
        summary_lines = []
        length_threshold_lines = []
        quality_threshold_lines = []

        for line in f["f"]:
            line = line.strip()
            if line.startswith("Nanoq Read Summary"):
                segment = "summary"
                continue
            elif line.startswith("Read length thresholds"):
                segment = "length_thresholds"
                continue
            elif line.startswith("Read quality thresholds"):
                segment = "quality_thresholds"
                continue

            if segment == "summary":
                summary_lines.append(line)
            elif segment == "length_thresholds":
                length_threshold_lines.append(line)
            elif segment == "quality_thresholds":
                quality_threshold_lines.append(line)

        # Parse each part
        nanoq_summary = parse_summary(summary_lines)
        read_length_thresholds = parse_thresholds(length_threshold_lines, "bp")
        read_quality_thresholds = parse_thresholds(quality_threshold_lines, "")

        # Prepare the final nanoq stats dictionary
        nanoq_stats: Dict[str, float] = {}
        for metric, value in zip(nanoq_summary["Metric"], nanoq_summary["Value"]):
            if metric in self._KEYS_MAPPING:
                nanoq_stats[self._KEYS_MAPPING[metric]] = value

        for threshold, num_reads in zip(read_length_thresholds["Threshold"], read_length_thresholds["Number of Reads"]):
            nanoq_stats[f"Reads > {threshold}"] = num_reads

        for threshold, num_reads in zip(
            read_quality_thresholds["Threshold"], read_quality_thresholds["Number of Reads"]
        ):
            nanoq_stats[f"Reads > Q{threshold}"] = num_reads

        # Store nanoq_summary for general stats
        self.nanoq_summary[f["s_name"]] = {k: v for k, v in zip(nanoq_summary["Metric"], nanoq_summary["Value"])}

        self.save_data(f, nanoq_stats)

    def save_data(self, f, nanoq_stats: Dict[str, float]) -> None:
        """
        Normalise fields and save parsed data.
        """
        if any(k.startswith("Reads > Q") for k in nanoq_stats):
            stat_type = "quality"
        else:
            stat_type = "summary"

        out_d = {f"{k}_{stat_type}": v for k, v in nanoq_stats.items()}

        # Warn if we find overlapping data for the same sample
        if f["s_name"] in self.nanoq_data:
            if not set(self.nanoq_data[f["s_name"]].keys()).isdisjoint(out_d.keys()):
                log.debug(f"Duplicate sample data found! Overwriting: {f['s_name']}")

        self.nanoq_data.setdefault(f["s_name"], {}).update(out_d)

        self.add_data_source(f)

    def nanoq_general_stats(self) -> None:
        """Nanoq General Stats Table"""
        headers = {
            "Number of reads": {
                "title": f"# Reads ({config.long_read_count_prefix})",
                "description": f"Number of reads ({config.long_read_count_desc})",
                "scale": "YlGn",
                "shared_key": "long_read_count",
                "modify": lambda x: x * config.long_read_count_multiplier,
            },
            "Number of bases": {
                "title": f"Total Bases ({config.base_count_prefix})",
                "description": f"Total bases ({config.base_count_desc})",
                "scale": "BuPu",
                "shared_key": "base_count",
                "modify": lambda x: x * config.base_count_multiplier,
            },
            "N50 read length": {
                "title": "Read N50",
                "description": "N50 read length",
                "scale": "RdPu",
                "suffix": " bp",
            },
            "Longest read": {
                "title": "Longest Read",
                "description": "Longest read length (bp)",
                "scale": "PuBuGn",
                "suffix": " bp",
                "shared_key": "nucleotides",
            },
            "Shortest read": {
                "title": "Shortest Read",
                "description": "Shortest read length (bp)",
                "scale": "PuBuGn",
                "suffix": " bp",
                "shared_key": "nucleotides",
            },
            "Mean read length": {
                "title": "Mean Length",
                "description": "Mean read length (bp)",
                "scale": "Purples",
                "suffix": " bp",
                "shared_key": "nucleotides",
            },
            "Median read length": {
                "title": "Median Length",
                "description": "Median read length (bp)",
                "scale": "BuGn",
                "suffix": " bp",
                "shared_key": "nucleotides",
            },
            "Mean read quality": {
                "title": "Mean Qual",
                "description": "Mean read quality (Phred scale)",
                "scale": "PiYG",
                "shared_key": "phred_score",
            },
            "Median read quality": {
                "title": "Median Qual",
                "description": "Median read quality (Phred scale)",
                "scale": "RdYlGn",
                "shared_key": "phred_score",
            },
        }

        # Add columns to the general stats table
        self.general_stats_addcols(self.nanoq_summary, headers)

    def reads_by_quality_plot(self) -> None:
        def _get_total_reads(_data_dict: Dict[str, Any]) -> Any:
            for _stat_type in self._stat_types:
                total_key = f"Number of reads_{_stat_type}"
                if total_key in _data_dict:
                    return _data_dict[total_key]
            return None

        # Get data for plot
        data: Dict[str, Dict[str, float]] = {}
        cats: List[str] = []
        for name, d in self.nanoq_data.items():
            # Iterate over available types and keep categories in same plot
            for stat_type in self._stat_types:
                reads_by_q = {k: v for k, v in d.items() if re.match(r"Reads > Q\d+_" + stat_type, k)}
                if len(reads_by_q) == 0:
                    continue

                data.setdefault(name, {})

                total_reads = _get_total_reads(d)
                for k, v in reads_by_q.items():
                    threshold = k.split("_")[0]
                    data[name][threshold] = v / total_reads * 100
                    cats.append(threshold)

        if len(data) == 0:
            return None

        # Plot
        cats = sorted(cats, key=lambda x: int(re.sub(r"[^\d.]", "", x)))
        pconfig = {
            "id": "nanoq_plot_quality",
            "title": "Read Quality",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }
        self.add_section(
            name="Reads Quality",
            anchor="nanoq-plot-quality",
            description="Read counts categorised by read quality (Phred score).",
            helptext="""
            Sequencing machines assign each generated read a quality score using the
            [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score).
            The phred score represents the liklelyhood that a given read contains errors.
            High quality reads have a high score.
            """,
            plot=bargraph.plot(data, cats, pconfig),
        )

    def reads_by_length_plot(self) -> None:
        def _get_total_reads(_data_dict: Dict[str, Any]) -> Any:
            for _stat_type in self._stat_types:
                total_key = f"Number of reads_{_stat_type}"
                if total_key in _data_dict:
                    return _data_dict[total_key]
            return None

        # Get data for plot
        data: Dict[str, Dict[str, float]] = {}
        cats: List[str] = []
        for name, d in self.nanoq_data.items():
            # Iterate over available types and keep categories in same plot
            for stat_type in self._stat_types:
                reads_by_l = {k: v for k, v in d.items() if re.match(r"Reads > \d+bp_" + stat_type, k)}
                if len(reads_by_l) == 0:
                    continue

                data.setdefault(name, {})

                total_reads = _get_total_reads(d)
                for k, v in reads_by_l.items():
                    threshold = k.split("_")[0]
                    data[name][threshold] = v / total_reads * 100
                    cats.append(threshold)

        if len(data) == 0:
            return None

        # Plot
        cats = sorted(cats, key=lambda x: int(re.sub(r"[^\d.]", "", x)))
        pconfig = {
            "id": "nanoq_plot_length",
            "title": "Read Length",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(
            name="Reads Length",
            anchor="nanoq-plot-length",
            description="Read counts categorised by read length.",
            plot=bargraph.plot(data, cats, pconfig),
        )
