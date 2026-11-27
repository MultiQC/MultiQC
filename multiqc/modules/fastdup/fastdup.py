import logging
import math
import re
from typing import Dict, Optional, Tuple

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.types import SectionAlert

log = logging.getLogger(__name__)

COUNT_FIELDS = (
    "UNPAIRED_READS_EXAMINED",
    "READ_PAIRS_EXAMINED",
    "SECONDARY_OR_SUPPLEMENTARY_RDS",
    "UNMAPPED_READS",
    "UNPAIRED_READ_DUPLICATES",
    "READ_PAIR_DUPLICATES",
    "READ_PAIR_OPTICAL_DUPLICATES",
    "ESTIMATED_LIBRARY_SIZE",
)


def parse_report(contents: str) -> Tuple[Dict, Dict[int, Optional[float]]]:
    """Parse the single aggregate metrics row and predicted-yield histogram."""
    metrics: Dict = {}
    histogram: Dict[int, Optional[float]] = {}
    keys = []
    block = None
    for line in contents.splitlines():
        if line.strip() == "## METRICS":
            if metrics:
                raise ValueError("Multiple FastDup metrics tables in one report")
            block = "metrics_header"
        elif line.startswith("## HISTOGRAM"):
            block = "histogram_header"
        elif not line.strip():
            block = None
        elif line.startswith("#"):
            continue
        elif block == "metrics_header":
            keys = line.split("\t")
            required = {"LIBRARY", "PERCENT_DUPLICATION", *COUNT_FIELDS}
            if not required.issubset(keys) or len(keys) != len(set(keys)):
                raise ValueError("Missing or duplicate FastDup metrics columns")
            block = "metrics"
        elif block == "metrics":
            if metrics:
                raise ValueError("Expected one aggregate FastDup metrics row")
            values = line.split("\t")
            if len(values) != len(keys):
                raise ValueError("FastDup metrics row does not match its header")
            metrics = dict(zip(keys, values))
            for field in COUNT_FIELDS:
                metrics[field] = int(metrics[field])
                if metrics[field] < 0:
                    raise ValueError(f"Negative FastDup count: {field}")
            metrics["PERCENT_DUPLICATION"] = float(metrics["PERCENT_DUPLICATION"])
            if not math.isfinite(metrics["PERCENT_DUPLICATION"]) or not 0 <= metrics["PERCENT_DUPLICATION"] <= 1:
                raise ValueError("Invalid FastDup duplication fraction")
            if (
                metrics["READ_PAIR_OPTICAL_DUPLICATES"] > metrics["READ_PAIR_DUPLICATES"]
                or metrics["READ_PAIR_DUPLICATES"] > metrics["READ_PAIRS_EXAMINED"]
                or metrics["UNPAIRED_READ_DUPLICATES"] > metrics["UNPAIRED_READS_EXAMINED"]
            ):
                raise ValueError("FastDup duplicate counts exceed examined counts")
        elif block == "histogram_header":
            # FastDup writes a space here, but tabs in the data rows.
            if line.split() != ["BIN", "CoverageMult"]:
                raise ValueError("Unexpected FastDup histogram columns")
            block = "histogram"
        elif block == "histogram":
            fields = line.split()
            if len(fields) != 2:
                raise ValueError("Invalid FastDup histogram row")
            depth, value = int(fields[0]), float(fields[1])
            if depth <= 0 or depth in histogram:
                raise ValueError("Invalid or duplicate FastDup histogram bin")
            if math.isfinite(value) and value < 0:
                raise ValueError("Negative FastDup predicted yield")
            # FastDup can emit NaN when library size cannot be estimated.
            histogram[depth] = value if math.isfinite(value) else None

    if not metrics:
        raise ValueError("Missing FastDup metrics row")
    if block == "histogram_header":
        raise ValueError("Missing FastDup histogram header")
    return metrics, histogram


class MultiqcModule(BaseMultiqcModule):
    """
    Parses the metrics file produced by `fastdup --metrics sample.metrics.txt`.
    The duplication fraction, duplicate-pair counts and optical duplicate-pair
    counts are available in General Statistics. Pair-count columns are hidden
    by default and can be enabled with the table's Configure Columns control.

    Sample names come from `--input` in the recorded command line, falling back
    to the metrics filename. To always use filenames, set
    `use_filename_as_sample_name: [fastdup]`.

    The histogram describes predicted unique-read yield at increased sequencing
    depth, rather than observed duplicate frequencies. A curve is unavailable
    when the library size cannot be estimated.
    """

    def __init__(self):
        super().__init__(
            name="FastDup",
            anchor="fastdup",
            href="https://github.com/zzhofict/FastDup",
            info="Identifies duplicate reads in alignment files.",
            doi="10.1093/bioinformatics/btaf633",
            license="MIT License",
            license_url="https://github.com/zzhofict/FastDup/blob/main/LICENSE",
        )
        self.metrics: Dict[str, Dict] = {}
        self.histograms: Dict[str, Dict[int, Optional[float]]] = {}

        for f in self.find_log_files("fastdup"):
            try:
                metrics, histogram = parse_report(f["f"])
            except ValueError as error:
                raise ValueError(f"Invalid FastDup report {f['root']}/{f['fn']}: {error}") from error

            s_name = f["s_name"]
            # Require the path to end at another option or the end of the line.
            # This avoids taking just the first word of an unquoted path with spaces.
            match = re.search(
                r""" --input(?:=|\s+)("[^"]+"|'[^']+'|\S+)(?=\s+-\S|\s*$)""",
                f["f"],
                flags=re.MULTILINE,
            )
            if match:
                input_name = re.split(r"[/\\]", match.group(1).strip("\"'"))[-1]
                s_name = self.clean_s_name(input_name, f)

            if s_name in self.metrics:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.metrics[s_name] = metrics
            self.histograms[s_name] = histogram
            self.add_data_source(f, s_name)
            self.add_software_version(None, s_name)

        self.metrics = self.ignore_samples(self.metrics)
        self.histograms = {s: self.histograms[s] for s in self.metrics}
        if not self.metrics:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.metrics)} reports")

        self.general_stats_addcols(
            self.metrics,
            {
                "PERCENT_DUPLICATION": {
                    "title": "Duplication",
                    "description": "Fraction of examined reads marked as duplicates",
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "scale": "OrRd",
                    "modify": lambda value: value * 100,
                },
                "READ_PAIR_DUPLICATES": {
                    "title": "Duplicate Pairs",
                    "description": "Duplicate read pairs, including optical duplicates",
                    "scale": "Purples",
                    "hidden": True,
                    "format": "{:,.0f}",
                },
                "READ_PAIR_OPTICAL_DUPLICATES": {
                    "title": "Optical Duplicate Pairs",
                    "description": "Optical duplicate read pairs, a subset of duplicate pairs",
                    "scale": "Oranges",
                    "hidden": True,
                    "format": "{:,.0f}",
                },
            },
        )
        self._add_duplication_section()
        self._add_yield_section()
        self.write_data_file(self.metrics, "multiqc_fastdup")
        self.write_data_file(self.histograms, "multiqc_fastdup_histogram")

    def _add_duplication_section(self):
        read_counts = {
            sample: {
                "unique_pairs": 2 * (m["READ_PAIRS_EXAMINED"] - m["READ_PAIR_DUPLICATES"]),
                "unique_unpaired": m["UNPAIRED_READS_EXAMINED"] - m["UNPAIRED_READ_DUPLICATES"],
                "optical_duplicates": 2 * m["READ_PAIR_OPTICAL_DUPLICATES"],
                "other_pair_duplicates": 2 * (m["READ_PAIR_DUPLICATES"] - m["READ_PAIR_OPTICAL_DUPLICATES"]),
                "unpaired_duplicates": m["UNPAIRED_READ_DUPLICATES"],
                "unmapped": m["UNMAPPED_READS"],
            }
            for sample, m in self.metrics.items()
        }
        plot_data = {sample: counts for sample, counts in read_counts.items() if any(counts.values())}
        empty_samples = [sample for sample in read_counts if sample not in plot_data]
        self.add_section(
            name="Duplication",
            anchor=f"{self.anchor}-duplication",
            description="Reads grouped by duplication state. Pair counts are doubled to show individual reads.",
            helptext="Optical duplicate pairs are a subset of duplicate pairs and are shown separately. "
            "Secondary and supplementary alignments are excluded from this plot. "
            "Raw pair counts are preserved in the exported metrics.",
            alerts=[SectionAlert(message="No reads available to plot.", affected_samples=empty_samples)]
            if empty_samples
            else [],
            plot=bargraph.plot(
                plot_data,
                {
                    "unique_pairs": {"name": "Unique paired reads"},
                    "unique_unpaired": {"name": "Unique unpaired reads"},
                    "optical_duplicates": {"name": "Optical duplicate paired reads"},
                    "other_pair_duplicates": {"name": "Other duplicate paired reads"},
                    "unpaired_duplicates": {"name": "Duplicate unpaired reads"},
                    "unmapped": {"name": "Unmapped reads"},
                },
                {
                    "id": f"{self.anchor}_duplication",
                    "title": "FastDup: Duplication",
                    "ylab": "Reads",
                    "cpswitch_c_active": False,
                },
            )
            if plot_data
            else None,
        )

    def _add_yield_section(self):
        plot_data = {}
        unavailable_samples = []
        nonfinite_samples = []
        for sample, histogram in self.histograms.items():
            metrics = self.metrics[sample]
            finite_values = {depth: value for depth, value in histogram.items() if value is not None}
            if len(finite_values) != len(histogram):
                nonfinite_samples.append(sample)
            if (
                finite_values
                and metrics["ESTIMATED_LIBRARY_SIZE"] > 0
                and metrics["READ_PAIRS_EXAMINED"] > metrics["READ_PAIR_DUPLICATES"]
            ):
                plot_data[sample] = finite_values
            else:
                unavailable_samples.append(sample)
        alerts = []
        if unavailable_samples:
            alerts.append(
                SectionAlert(
                    message="Predicted yield is unavailable because the histogram is absent or "
                    "library size could not be estimated.",
                    affected_samples=unavailable_samples,
                )
            )
        if nonfinite_samples:
            alerts.append(
                SectionAlert(
                    message="Nonfinite predicted yield values were excluded from the plot.",
                    affected_samples=nonfinite_samples,
                )
            )
        self.add_section(
            name="Predicted Yield",
            anchor=f"{self.anchor}-yield",
            description="Predicted increase in unique-read yield when sequencing the same library more deeply.",
            helptext="The curve uses the CoverageMult values reported by FastDup, without recalculation. "
            "The x-axis is sequencing depth relative to the current depth; the y-axis is predicted "
            "unique-read yield relative to the currently observed unique read pairs.",
            alerts=alerts,
            plot=linegraph.plot(
                plot_data,
                {
                    "id": f"{self.anchor}_yield",
                    "title": "FastDup: Predicted Unique-Read Yield",
                    "xlab": "Sequencing-depth multiple",
                    "ylab": "Predicted unique-yield multiple",
                    "xmin": 1,
                    "ymin": 0,
                },
            )
            if plot_data
            else None,
        )
