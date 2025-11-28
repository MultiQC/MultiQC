"""MultiQC submodule to parse output from seqkit stats"""

import logging
from typing import Dict, Optional

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


def parse_seqkit_stats(module: BaseMultiqcModule) -> int:
    """Find seqkit stats logs and parse their data"""

    seqkit_stats: Dict[str, Dict] = {}

    for f in module.find_log_files("seqkit/stats"):
        parsed_data = parse_stats_report(f["f"], f["s_name"])
        if parsed_data:
            for s_name, data in parsed_data.items():
                s_name = module.clean_s_name(s_name, f)
                if s_name in seqkit_stats:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                module.add_data_source(f, s_name=s_name, section="stats")
                seqkit_stats[s_name] = data

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                module.add_software_version(None, s_name)

    # Filter to strip out ignored sample names
    seqkit_stats = module.ignore_samples(seqkit_stats)

    if len(seqkit_stats) == 0:
        return 0

    # General Stats Table - MultiQC automatically filters missing columns
    general_stats_headers: Dict = {
        "num_seqs": {
            "title": "# Seqs",
            "description": "Number of sequences",
            "scale": "Blues",
            "format": "{:,.0f}",
            "shared_key": "seq_count",
        },
        "sum_len": {
            "title": "Total bp",
            "description": "Total number of bases",
            "scale": "Greens",
            "format": "{:,.0f}",
            "shared_key": "base_count",
            "hidden": True,
        },
        "avg_len": {
            "title": "Avg Length",
            "description": "Average sequence length",
            "scale": "Purples",
            "format": "{:,.1f}",
            "suffix": " bp",
        },
        "GC_pct": {
            "title": "GC%",
            "description": "GC content percentage",
            "min": 0,
            "max": 100,
            "scale": "Greens",
            "suffix": "%",
            "hidden": True,
        },
        "Q30_pct": {
            "title": "Q30%",
            "description": "Percentage of bases with quality score >= 30",
            "min": 0,
            "max": 100,
            "scale": "RdYlGn",
            "suffix": "%",
        },
        "AvgQual": {
            "title": "Avg Qual",
            "description": "Average quality score",
            "scale": "RdYlGn",
            "hidden": True,
        },
        "N50": {
            "title": "N50",
            "description": "N50 sequence length",
            "scale": "Oranges",
            "format": "{:,.0f}",
            "suffix": " bp",
            "hidden": True,
        },
    }

    # Get general stats headers using the utility function
    stats_headers = module.get_general_stats_headers(all_headers=general_stats_headers)

    # Add headers to general stats table
    if stats_headers:
        module.general_stats_addcols(seqkit_stats, stats_headers, namespace="seqkit")

    # Create detailed table with all columns
    table_headers: Dict = {
        "file": {
            "title": "File",
            "description": "Input file name",
            "hidden": True,
        },
        "format": {
            "title": "Format",
            "description": "File format (FASTA or FASTQ)",
        },
        "type": {
            "title": "Type",
            "description": "Sequence type (DNA, RNA, Protein)",
        },
        "num_seqs": {
            "title": "# Seqs",
            "description": "Number of sequences",
            "format": "{:,.0f}",
            "shared_key": "seq_count",
        },
        "sum_len": {
            "title": "Total Length",
            "description": "Total number of bases/residues",
            "format": "{:,.0f}",
            "shared_key": "base_count",
        },
        "min_len": {
            "title": "Min Length",
            "description": "Minimum sequence length",
            "format": "{:,.0f}",
        },
        "avg_len": {
            "title": "Avg Length",
            "description": "Average sequence length",
            "format": "{:,.1f}",
        },
        "max_len": {
            "title": "Max Length",
            "description": "Maximum sequence length",
            "format": "{:,.0f}",
        },
        "Q1": {
            "title": "Q1",
            "description": "First quartile of sequence length",
            "format": "{:,.0f}",
        },
        "Q2": {
            "title": "Q2 (Median)",
            "description": "Median sequence length",
            "format": "{:,.0f}",
        },
        "Q3": {
            "title": "Q3",
            "description": "Third quartile of sequence length",
            "format": "{:,.0f}",
        },
        "sum_gap": {
            "title": "Gaps",
            "description": "Total number of gaps",
            "format": "{:,.0f}",
        },
        "N50": {
            "title": "N50",
            "description": "N50 sequence length",
            "format": "{:,.0f}",
        },
        "N50_num": {
            "title": "N50 Count",
            "description": "Number of sequences >= N50",
            "format": "{:,.0f}",
        },
        "Q20_pct": {
            "title": "Q20%",
            "description": "Percentage of bases with quality score >= 20",
            "min": 0,
            "max": 100,
            "suffix": "%",
        },
        "Q30_pct": {
            "title": "Q30%",
            "description": "Percentage of bases with quality score >= 30",
            "min": 0,
            "max": 100,
            "suffix": "%",
        },
        "AvgQual": {
            "title": "Avg Quality",
            "description": "Average quality score",
        },
        "GC_pct": {
            "title": "GC%",
            "description": "GC content percentage",
            "min": 0,
            "max": 100,
            "suffix": "%",
        },
        "sum_n": {
            "title": "N Bases",
            "description": "Total number of N bases",
            "format": "{:,.0f}",
        },
    }

    # Add table section
    module.add_section(
        name="Stats",
        anchor="seqkit-stats",
        description="Statistics from <code>seqkit stats</code> showing sequence file metrics.",
        plot=table.plot(
            seqkit_stats,
            table_headers,
            pconfig={
                "id": "seqkit-stats-table",
                "title": "SeqKit: Stats",
                "namespace": "seqkit",
            },
        ),
    )

    # Create bar plot for sequence counts and lengths
    bargraph_data_seqs = {s_name: {"Sequences": d.get("num_seqs", 0)} for s_name, d in seqkit_stats.items()}

    bargraph_data_bases = {s_name: {"Total Bases": d.get("sum_len", 0)} for s_name, d in seqkit_stats.items()}

    module.add_section(
        name="Sequence Counts",
        anchor="seqkit-stats-seqcounts",
        description="Number of sequences per sample.",
        plot=bargraph.plot(
            bargraph_data_seqs,
            pconfig={
                "id": "seqkit-stats-seq-counts",
                "title": "SeqKit: Sequence Counts",
                "ylab": "Number of Sequences",
            },
        ),
    )

    module.add_section(
        name="Total Bases",
        anchor="seqkit-stats-bases",
        description="Total number of bases per sample.",
        plot=bargraph.plot(
            bargraph_data_bases,
            pconfig={
                "id": "seqkit-stats-total-bases",
                "title": "SeqKit: Total Bases",
                "ylab": "Number of Bases",
            },
        ),
    )

    # Write parsed report data to a file
    module.write_data_file(seqkit_stats, "multiqc_seqkit_stats")

    return len(seqkit_stats)


def parse_stats_report(file_content: str, fallback_sample_name: Optional[str] = None) -> Dict[str, Dict]:
    """
    Parse seqkit stats output file.

    Returns a dictionary with sample names as keys and parsed data as values.
    Handles both tab-separated (--tabular flag) and space-separated (default) output.
    """
    parsed_data: Dict[str, Dict] = {}
    lines = file_content.strip().split("\n")

    if len(lines) < 2:
        return {}

    # Determine delimiter: tab-separated or space-separated
    header_line = lines[0]
    delimiter = "\t" if "\t" in header_line else None
    headers = header_line.split(delimiter)
    headers = [h.strip() for h in headers]

    # Check if this looks like seqkit stats output
    expected_headers = ["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"]
    if not all(h in headers for h in expected_headers):
        return {}

    # Define column types for parsing
    int_columns = {"num_seqs", "sum_len", "min_len", "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "sum_n"}
    float_columns = {"avg_len", "Q20(%)", "Q30(%)", "AvgQual", "GC(%)"}
    # Mapping for columns with special characters in names
    column_renames = {"Q20(%)": "Q20_pct", "Q30(%)": "Q30_pct", "GC(%)": "GC_pct"}

    # Parse data lines
    for line in lines[1:]:
        if not line.strip():
            continue

        values = line.split(delimiter)
        if len(values) < len(headers):
            continue

        row = dict(zip(headers, [v.strip() for v in values]))
        data: Dict = {}

        # Get sample name from file column
        file_value = row.get("file", "")
        if file_value == "-":
            # Use fallback sample name for stdin input
            sample_name = fallback_sample_name if fallback_sample_name else "stdin"
        else:
            # Handle both Unix and Windows path separators
            sample_name = file_value.replace("\\", "/").split("/")[-1]
            # Remove common sequence file extensions
            for ext in [".fq.gz", ".fastq.gz", ".fq", ".fastq", ".fa.gz", ".fasta.gz", ".fa", ".fasta"]:
                if sample_name.endswith(ext):
                    sample_name = sample_name[: -len(ext)]
                    break

        data["file"] = file_value

        # Parse remaining columns
        for header, value in row.items():
            if header == "file":
                continue
            elif header in ("format", "type"):
                data[header] = value
            elif header in int_columns:
                try:
                    data[header] = int(value.replace(",", ""))
                except (ValueError, AttributeError):
                    pass
            elif header in float_columns:
                try:
                    key = column_renames.get(header, header)
                    data[key] = float(value.replace(",", ""))
                except (ValueError, AttributeError):
                    pass

        if sample_name and data:
            parsed_data[sample_name] = data

    return parsed_data
