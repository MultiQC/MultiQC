import logging
import re
from collections import defaultdict
from copy import deepcopy
from typing import Dict, List, Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table, linegraph
from multiqc.plots.table_object import TableConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    _stat_types = ("summary", "quality")

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="nanoq",
            anchor="nanoq",
            href="https://github.com/nerdna/nanoq/",
            info="Reports read quality and length from nanopore sequencing data",
            doi="10.21105/joss.02991",
        )

        # Find and load any nanoq reports
        data_by_sample: Dict[str, Dict[str, float]] = {}
        for f in self.find_log_files("nanoq", filehandles=True):
            sample_data = parse_nanoq_log(f)
            if sample_data:
                data_by_sample[f["s_name"]] = sample_data
                if f["s_name"] in data_by_sample:
                    log.debug(f"Duplicate sample data found! Overwriting: {f['s_name']}")

                self.add_data_source(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        # Write parsed report data to a file
        self.write_data_file(data_by_sample, "multiqc_nanoq")

        # Add Nanoq summary to the general stats table
        self.add_table(data_by_sample)

        # Quality distribution Plot
        self.reads_by_quality_plot(data_by_sample)

        # Read length distribution Plot
        self.reads_by_length_plot(data_by_sample)

    def add_table(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        headers: Dict[str, Dict] = {
            "Number of reads": {
                "title": "Reads",
                "description": "Number of reads",
                "scale": "Greens",
                "shared_key": "read_count",
            },
            "Number of bases": {
                "title": "Bases",
                "description": "Total bases sequenced",
                "scale": "Purples",
                "shared_key": "base_count",
            },
            "N50 read length": {
                "title": "Read N50",
                "description": "N50 read length",
                "scale": "Blues",
                "suffix": "bp",
                "format": "{:,.0f}",
            },
            "Longest read": {
                "title": "Longest Read",
                "description": "Longest read length",
                "suffix": "bp",
                "scale": "Oranges",
                "format": "{:,.0f}",
            },
            "Shortest read": {
                "title": "Shortest Read",
                "description": "Shortest read length",
                "suffix": "bp",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
            "Mean read length": {
                "title": "Mean Length",
                "description": "Mean read length",
                "suffix": "bp",
                "scale": "PuBuGn",
            },
            "Median read length": {
                "title": "Median Length",
                "description": "Median read length (bp)",
                "scale": "RdYlBu",
                "format": "{:,.0f}",
            },
            "Mean read quality": {
                "title": "Mean Qual",
                "description": "Mean read quality (Phred scale)",
                "scale": "PiYG",
            },
            "Median read quality": {
                "title": "Median Qual",
                "description": "Median read quality (Phred scale)",
                "scale": "Spectral",
            },
        }

        self.add_section(
            name="Nanoq Summary",
            anchor="nanoq-summary",
            description="Statistics from Nanoq reports",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig=TableConfig(
                    id="nanoq_table",
                    title="Nanoq Summary",
                ),
            ),
        )

        general_stats_headers = deepcopy(headers)
        for k, h in general_stats_headers.items():
            h["hidden"] = True
        general_stats_headers["Number of reads"]["hidden"] = False
        general_stats_headers["N50 read length"]["hidden"] = False

        # Add columns to the general stats table
        self.general_stats_addcols(data_by_sample, general_stats_headers)

    def reads_by_quality_plot(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        # Get data for plot
        lineplot_data: Dict[str, Dict[int, float]] = defaultdict(dict)
        for name, d in data_by_sample.items():
            reads_by_q = {k: v for k, v in d.items() if k.startswith("Reads > Q")}
            if len(reads_by_q) == 0:
                continue

            total_reads = d["Number of reads"]
            for k, v in reads_by_q.items():
                threshold = int(k.split("> Q")[1])
                lineplot_data[name][threshold] = v / total_reads * 100

        # Plot
        self.add_section(
            name="Read quality",
            anchor="nanoq_plot_quality",
            description="Read counts categorised by read quality (Phred score).",
            helptext="""
            Sequencing machines assign each generated read a quality score using the
            [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score).
            The phred score represents the liklelyhood that a given read contains errors.
            High quality reads have a high score.
            """,
            plot=linegraph.plot(
                lineplot_data,
                linegraph.LinePlotConfig(
                    id="nanoq_plot_quality_plot",
                    title="Nanoq: read qualities",
                ),
            ),
        )

    def reads_by_length_plot(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        # Get data for plot
        linegraph_data: Dict[str, Dict[int, float]] = defaultdict(dict)
        for name, d in data_by_sample.items():
            reads_by_l = {k: v for k, v in d.items() if re.match(r"Reads > \d+bp", k)}
            if len(reads_by_l) == 0:
                continue

            for k, v in reads_by_l.items():
                threshold = int(k.split(" > ")[1].strip("bp"))
                linegraph_data[name][threshold] = v

        self.add_section(
            name="Read lengths",
            anchor="nanoq_plot_length",
            description="Read counts categorised by read length.",
            plot=linegraph.plot(
                linegraph_data,
                linegraph.LinePlotConfig(
                    id="nanoq_plot_length_plot",
                    title="Nanoq: read lengths",
                    xlog=True,
                ),
            ),
        )


def parse_nanoq_log(f) -> Dict[str, float]:
    """Parse output from nanoq"""
    stats: Dict[str, float] = dict()

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

    for line in summary_lines:
        if ":" in line:
            metric, value = line.split(":", 1)
            stats[metric.strip()] = float(value.strip())

    # Helper function to parse thresholds part
    def parse_thresholds(lines: List[str], threshold_type: str) -> Dict[str, List[Any]]:
        _thresholds: Dict[str, List[Any]] = {"Threshold": [], "Number of Reads": [], "Percentage": []}
        for _line in lines:
            match = re.match(r">\s*(\d+)\s+(\d+)\s+(\d+\.\d+)%", _line)
            if match:
                _threshold, _num_reads, _percentage = match.groups()
                _thresholds["Threshold"].append(_threshold + (threshold_type if threshold_type == "bp" else ""))
                _thresholds["Number of Reads"].append(int(_num_reads))
                _thresholds["Percentage"].append(float(_percentage))
        return _thresholds

    read_length_thresholds = parse_thresholds(length_threshold_lines, "bp")
    read_quality_thresholds = parse_thresholds(quality_threshold_lines, "")

    for threshold, num_reads in zip(read_length_thresholds["Threshold"], read_length_thresholds["Number of Reads"]):
        stats[f"Reads > {threshold}"] = num_reads

    for threshold, num_reads in zip(read_quality_thresholds["Threshold"], read_quality_thresholds["Number of Reads"]):
        stats[f"Reads > Q{threshold}"] = num_reads

    return stats
