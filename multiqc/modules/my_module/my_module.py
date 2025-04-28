import logging
import re
from collections import defaultdict
from copy import deepcopy
from typing import Callable, Dict, List, Any, Tuple, Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table, bargraph
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table_object import TableConfig, ColumnMeta
from multiqc.utils import mqc_colour
from multiqc import config

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Qualalyser provides multiple subcommands, and the MultiQC module currently only supports `quality`.

    Qualalyser outputs useful information into stdout, and you need to capture it to
    a file for the module to recognize. To pipe stderr into a file, run the tool
    as follows:

    qualalyser quality 2> sample1.log

    Note the that the sample name is parsed from the filename by default, in this case,
    the reported name will be "sample1".

    #### Configuration

    By default, Qualalyser uses the following quality threshold: 10.

    To override it, use the following config:

    qualalyser:
      min_quality: 10

    Version 1.1.0 of Qualalyser is tested.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Qualalyser",
            anchor="qualalyser",
            href="https://github.com/bioinformatics-centre/qualalyser/",
            info="Reports read quality and length from sequencing data",
            doi="10.21105/joss.02991",
        )

        # Find and load any Qualalyser reports
        data_by_sample: Dict[str, Dict[str, float]] = {}
        for f in self.find_log_files("qualalyser/quality", filehandles=True):
            sample_data = parse_qualalyser_log(f)
            if sample_data:
                s_name = f['s_name']
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                data_by_sample[s_name] = sample_data
                self.add_data_source(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Add Qualalyser summary to the general stats table
        self.add_table(data_by_sample)

        # Quality distribution Plot
        self.reads_by_quality_plot(data_by_sample)

        # Read length distribution Plot
        self.reads_by_length_plot(data_by_sample)

        # Write parsed report data to a file
        self.write_data_file(data_by_sample, "multiqc_qualalyser")

    def add_table(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        headers: Dict[str, Dict] = {
            "Number of reads": ColumnMeta(
                title="Reads",
                description="Number of reads",
                scale="Greens",
                shared_key="read_count",
            ),
            "Number of bases": ColumnMeta(
                title="Bases",
                description="Total bases sequenced",
                scale="Purples",
                shared_key="base_count",
            ),
            "N50 read length": ColumnMeta(
                title="Read N50",
                description="N50 read length",
                scale="Blues",
                suffix="bp",
                format="{:,.0f}",
            ),
            "Longest read": ColumnMeta(
                title="Longest Read",
                description="Longest read length",
                suffix="bp",
                scale="Oranges",
                format="{:,.0f}",
            ),
            "Mean read length": ColumnMeta(
                title="Mean Length",
                description="Mean read length",
                suffix="bp",
                scale="PuBuGn",
            ),
            "Median read length": ColumnMeta(
                title="Median Length",
                description="Median read length (bp)",
                scale="RdYlBu",
                format="{:,.0f}",
            ),
            "Mean read quality": ColumnMeta(
                title="Mean Qual",
                description="Mean read quality (Phred scale)",
                scale="PiYG",
            ),
            "Median read quality": ColumnMeta(
                title="Median Qual",
                description="Median read quality (Phred scale)",
                scale="Spectral",
            ),
        }

        self.add_section(
            name="Qualalyser Summary",
            anchor="qualalyser-summary",
            description="Statistics from Qualalyser reports",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig=TableConfig(
                    id="qualalyser_table",
                    title="Qualalyser Summary",
                ),
            ),
        )

        # Add general stats table - hide all columns except for two
        general_stats_headers = deepcopy(headers)
        for h in general_stats_headers.values():
            h["hidden"] = True
        general_stats_headers["Number of reads"]["hidden"] = False
        general_stats_headers["N50 read length"]["hidden"] = False

        # Add columns to the general stats table
        self.general_stats_addcols(data_by_sample, general_stats_headers)

    def reads_by_quality_plot(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        barplot_data: Dict[str, Dict[str, float]] = defaultdict(dict)
        keys: List[str] = []
        min_quality = getattr(config, "qualalyser", {}).get("min_quality", 10)

        for name, d in data_by_sample.items():
            reads_by_q = {int(re.search(r"\d+", k).group(0)): v for k, v in d.items() if k.startswith("Reads > Q")}
            if not reads_by_q:
                continue

            thresholds = sorted(th for th in reads_by_q if th >= min_quality)
            if not thresholds:
                continue

            barplot_data[name], keys = get_ranges_from_cumsum(
                data=reads_by_q, thresholds=thresholds, total=d["Number of reads"], formatter=lambda x: f"Q{x}"
            )

        colours = mqc_colour.mqc_colour_scale("RdYlGn-rev", 0, len(keys))
        cats = {
            k: {"name": f"Reads {k}", "color": colours.get_colour(idx, lighten=1)} for idx, k in enumerate(keys[::-1])
        }

        # Plot
        self.add_section(
            name="Read quality",
            anchor="qualalyser_plot_quality",
            description="Read counts categorised by read quality (Phred score).",
            helptext="""
            Sequencing machines assign each generated read a quality score using the
            [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score).
            The phred score represents the liklelyhood that a given read contains errors.
            High quality reads have a high score.
            """,
            plot=bargraph.plot(
                barplot_data,
                cats,
                pconfig=bargraph.BarPlotConfig(
                    id="qualalyser_plot_quality_plot",
                    title="Qualalyser: read qualities",
                ),
            ),
        )


def parse_qualalyser_log(f) -> Dict[str, float]:
    """Parse output from Qualalyser"""
    stats: Dict[str, float] = dict()

    # Parse the file content
    segment = None
    summary_lines = []
    length_threshold_lines = []
    quality_threshold_lines = []

    for line in f["f"]:
        line = line.strip()
        if line.startswith("Qualalyser Read Summary"):
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

    return stats