"""MultiQC module to parse output from cramino"""

import fnmatch
import logging
from typing import Dict, Iterable, List, Tuple, Union, cast

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, table
from multiqc.plots.table_object import ColumnDict, ValueT

log = logging.getLogger(__name__)

KARYOTYPE_HEADER = "# Normalized read count per chromosome"
KARYOTYPE_EMPTY_WARNING = "# Warning - no contigs found in BAM file!"

# cramino prints "Median/Mean/Modal identity" for CIGAR-based identity, or "... est. identity"
# when an estimation mode was used instead; both represent the same metric.
IDENTITY_KEY_ALIASES = {
    "median_identity": ("Median identity", "Median est. identity"),
    "mean_identity": ("Mean identity", "Mean est. identity"),
    "modal_identity": ("Modal identity", "Modal est. identity"),
}

# Sensible default for the per-chromosome plot on a typical GRCh37/38 BAM: hide alt/decoy/
# random/unplaced contigs, which are noise for spotting whole-chromosome imbalances. Mirrors
# the example given in the mosdepth module's docstring for the same purpose.
DEFAULT_EXCLUDE_CONTIGS = ["*_alt", "*_decoy", "*_random", "chrUn*", "HLA*", "chrM", "chrEBV"]


def read_config() -> Dict:
    cfg = getattr(config, "cramino_config", dict())
    if not isinstance(cfg, dict):
        return {}

    cfg["include_contigs"] = cfg.get("include_contigs", [])
    if not isinstance(cfg["include_contigs"], list):
        cfg["include_contigs"] = []

    cfg["exclude_contigs"] = cfg.get("exclude_contigs", DEFAULT_EXCLUDE_CONTIGS)
    if not isinstance(cfg["exclude_contigs"], list):
        cfg["exclude_contigs"] = DEFAULT_EXCLUDE_CONTIGS

    return cfg


def parse_cramino_lines(lines: Iterable[str]) -> Tuple[Dict[str, str], Dict[str, float]]:
    """
    Parse cramino's default text output into (summary, karyotype).

    Format: `key<TAB>value` lines for file/alignment/read/identity/phase/splice stats, and an
    optional `# Normalized read count per chromosome` block of `chrom<TAB>normalized_count`
    lines (present only when cramino was run with `--karyotype`, and possibly followed by
    `# Warning - no contigs found in BAM file!` instead of data if none were found).
    """
    summary: Dict[str, str] = {}
    karyotype: Dict[str, float] = {}
    in_karyotype = False

    for raw_line in lines:
        stripped = raw_line.strip()
        if stripped == KARYOTYPE_HEADER:
            in_karyotype = True
            continue
        if stripped == KARYOTYPE_EMPTY_WARNING:
            in_karyotype = False
            continue
        if not stripped:
            # Blank lines are just block separators (including right after the karyotype
            # header itself); never carry data, so they can't be used to detect the end of
            # the karyotype block.
            continue

        fields = stripped.split("\t")
        if len(fields) != 2:
            raise ValueError(f"Expected key<TAB>value, got: {fields}")
        key, value = fields

        if in_karyotype:
            try:
                karyotype[key] = float(value)
                continue
            except ValueError:
                # Chromosome values are always numeric; a non-numeric value means we've
                # reached the trailing Path/Creation time pair, not more karyotype data.
                in_karyotype = False

        summary[key] = value

    return summary, karyotype


def build_cramino_stats(summary: Dict[str, str]) -> Dict[str, Union[int, float]]:
    """
    Extract the fields cramino always emits, plus optional identity stats (populated only for
    aligned input; absent when cramino was run on unaligned reads with `--ubam`).
    """
    stats: Dict[str, Union[int, float]] = {
        "num_alignments": int(summary["Number of alignments"]),
        "pct_from_total_alignments": float(summary["% from total alignments"]),
        "num_reads": int(summary["Number of reads"]),
        "yield_gb": float(summary["Yield [Gb]"]),
        "mean_coverage": float(summary["Mean coverage"]),
        "yield_gb_long": float(summary["Yield [Gb] (>25kb)"]),
        "n50": int(summary["N50"]),
        "n75": int(summary["N75"]),
        "median_length": float(summary["Median length"]),
        "mean_length": float(summary["Mean length"]),
    }

    for stat_key, aliases in IDENTITY_KEY_ALIASES.items():
        for alias in aliases:
            if alias in summary:
                stats[stat_key] = float(summary[alias])
                break

    return stats


def filter_karyotype_contigs(
    karyotype: Dict[str, float],
    include_contigs: List[str],
    exclude_contigs: List[str],
) -> Dict[str, float]:
    """Apply glob include/exclude patterns to a sample's per-chromosome data."""
    filtered: Dict[str, float] = {}
    for contig, value in karyotype.items():
        if exclude_contigs and any(fnmatch.fnmatch(contig, pattern) for pattern in exclude_contigs):
            continue
        if include_contigs and not any(fnmatch.fnmatch(contig, pattern) for pattern in include_contigs):
            continue
        filtered[contig] = value
    return filtered


class MultiqcModule(BaseMultiqcModule):
    """
    Cramino is a tool for quick quality assessment of aligned or unaligned long-read sequencing
    data (Oxford Nanopore or PacBio) in BAM or CRAM format.

    The module parses cramino's default text output (`cramino > sample.txt`), which always
    reports the number of reads, yield, read length (N50/N75, median, mean), and mean coverage.
    Identity statistics (median/mean/modal percent identity to the reference) are included when
    the input was aligned (not `--ubam`).

    When cramino was run with `--karyotype`, a normalized read count per chromosome is also
    included, and the module adds a per-chromosome plot useful for spotting chromosomal
    imbalances or off-target enrichment. By default, alt/decoy/random/unplaced contigs are
    hidden from that plot; this can be customised in the config file:

    ```yaml
    cramino_config:
      include_contigs:
        - "chr*"
      exclude_contigs:
        - "*_alt"
        - "*_decoy"
        - "*_random"
        - "chrUn*"
        - "HLA*"
        - "chrM"
        - "chrEBV"
    ```

    Note that exclusion takes precedence over inclusion.
    """

    def __init__(self):
        super().__init__(
            name="Cramino",
            anchor="cramino",
            href="https://github.com/wdecoster/cramino",
            info="Quick quality assessment of aligned or unaligned long-read sequencing data in BAM/CRAM format",
            doi="10.1093/bioinformatics/btad311",
        )

        self.cfg = read_config()

        cramino_data: Dict[str, Dict[str, Union[int, float]]] = {}
        karyotype_by_sample: Dict[str, Dict[str, float]] = {}

        for f in self.find_log_files("cramino", filehandles=True):
            s_name = self.clean_s_name(f["fn"], f)
            try:
                summary, karyotype = parse_cramino_lines(f["f"])
                stats = build_cramino_stats(summary)
            except (ValueError, KeyError) as e:
                log.warning(f"Could not parse cramino output in '{f['fn']}': {e}")
                continue

            if s_name in cramino_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name=s_name)
            cramino_data[s_name] = stats
            self.add_software_version(None, s_name)

            if karyotype:
                karyotype_by_sample[s_name] = karyotype

        cramino_data = self.ignore_samples(cramino_data)
        karyotype_by_sample = self.ignore_samples(karyotype_by_sample)

        if len(cramino_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(cramino_data)} reports")

        self.add_general_stats(cramino_data)
        self.add_summary_table(cramino_data)
        if karyotype_by_sample:
            self.add_karyotype_section(karyotype_by_sample)

        self.write_data_file(cramino_data, "multiqc_cramino")

    def add_general_stats(self, data: Dict[str, Dict[str, Union[int, float]]]) -> None:
        headers: Dict[str, ColumnDict] = {
            "yield_gb": {
                "title": "Yield (Gb)",
                "description": "Total sequencing yield in gigabases",
                "min": 0,
                "scale": "Blues",
            },
            "n50": {
                "title": "N50",
                "description": "Read length N50",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "RdPu",
            },
            "median_identity": {
                "title": "Median Identity",
                "description": "Median percent identity to the reference",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        }
        resolved_headers = self.get_general_stats_headers(all_headers=headers)
        if resolved_headers:
            self.general_stats_addcols(cast(Dict[str, Dict[str, ValueT]], data), resolved_headers)

    def add_summary_table(self, data: Dict[str, Dict[str, Union[int, float]]]) -> None:
        headers: Dict[str, ColumnDict] = {
            "num_reads": {
                "title": f"# Reads ({config.long_read_count_prefix})",
                "description": f"Number of reads ({config.long_read_count_desc}), excluding supplementary alignments",
                "modify": lambda x: x * config.long_read_count_multiplier,  # type: ignore[operator]
                "shared_key": "long_read_count",
                "scale": "YlGn",
            },
            "num_alignments": {
                "title": f"# Alignments ({config.long_read_count_prefix})",
                "description": f"Number of alignments used for the reported statistics ({config.long_read_count_desc})",
                "modify": lambda x: x * config.long_read_count_multiplier,  # type: ignore[operator]
                "shared_key": "long_read_count",
                "scale": "GnBu",
                "hidden": True,
            },
            "pct_from_total_alignments": {
                "title": "% Counted",
                "description": "Percentage of all alignment records in the file counted towards these statistics",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": True,
            },
            "yield_gb": {
                "title": "Yield (Gb)",
                "description": "Total sequencing yield in gigabases",
                "min": 0,
                "scale": "Blues",
            },
            "yield_gb_long": {
                "title": "Yield >25kb (Gb)",
                "description": "Sequencing yield from reads longer than 25kb, in gigabases",
                "min": 0,
                "scale": "PuBu",
                "hidden": True,
            },
            "mean_coverage": {
                "title": "Mean Cov.",
                "description": "Mean coverage (yield divided by genome size)",
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            },
            "n50": {
                "title": "N50",
                "description": "Read length N50",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "RdPu",
            },
            "n75": {
                "title": "N75",
                "description": "Read length N75",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "RdPu",
                "hidden": True,
            },
            "median_length": {
                "title": "Median Length",
                "description": "Median read length",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "BuPu",
                "hidden": True,
            },
            "mean_length": {
                "title": "Mean Length",
                "description": "Mean read length",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "Purples",
                "hidden": True,
            },
            "median_identity": {
                "title": "Median Identity",
                "description": "Median percent identity to the reference",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "mean_identity": {
                "title": "Mean Identity",
                "description": "Mean percent identity to the reference",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": True,
            },
            "modal_identity": {
                "title": "Modal Identity",
                "description": "Modal (most common) percent identity to the reference",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": True,
            },
        }

        self.add_section(
            name="Summary Statistics",
            anchor="cramino-stats",
            description="Summary statistics from <code>cramino</code>.",
            plot=table.plot(
                data,
                headers,
                pconfig={
                    "id": "cramino-stats-table",
                    "title": "Cramino: Summary Statistics",
                    "namespace": "Cramino",
                },
            ),
        )

    def add_karyotype_section(self, karyotype_by_sample: Dict[str, Dict[str, float]]) -> None:
        filtered_by_sample = {
            s_name: filter_karyotype_contigs(karyotype, self.cfg["include_contigs"], self.cfg["exclude_contigs"])
            for s_name, karyotype in karyotype_by_sample.items()
        }

        self.write_data_file(filtered_by_sample, "multiqc_cramino_karyotype")

        self.add_section(
            name="Normalized read count per chromosome",
            anchor="cramino-karyotype",
            description=(
                "Read count per chromosome, normalized to the sample's median, from <code>cramino --karyotype</code>."
            ),
            helptext="""
            A value close to 1 means the chromosome has roughly the same read depth as the rest
            of the genome. Values consistently above or below 1 across whole chromosomes can
            indicate aneuploidies (e.g. a 1.5x signal for a trisomy) or, for a targeted/enrichment
            experiment, off-target chromosomes.

            By default, alt/decoy/random/unplaced contigs are hidden; this can be customised with
            `cramino_config.include_contigs` / `cramino_config.exclude_contigs` in the MultiQC
            config file.
            """,
            plot=linegraph.plot(
                filtered_by_sample,
                pconfig={
                    "id": "cramino-karyotype-plot",
                    "title": "Cramino: Normalized read count per chromosome",
                    "xlab": "Chromosome",
                    "ylab": "Normalized read count",
                    "categories": True,
                },
            ),
        )
