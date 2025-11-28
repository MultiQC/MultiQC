"""MultiQC submodule to parse output from seqkit stats"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


def parse_seqkit_stats(module: BaseMultiqcModule) -> int:
    """Find seqkit stats logs and parse their data"""

    seqkit_stats: Dict[str, Dict] = {}

    for f in module.find_log_files("seqkit/stats"):
        parsed_data = parse_stats_report(f["f"])
        if parsed_data:
            for s_name, data in parsed_data.items():
                s_name = module.clean_s_name(s_name, f)
                if s_name in seqkit_stats:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                module.add_data_source(f, s_name=s_name, section="stats")
                seqkit_stats[s_name] = data

    # Filter to strip out ignored sample names
    seqkit_stats = module.ignore_samples(seqkit_stats)

    if len(seqkit_stats) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Check which columns are present in the data
    has_quality_metrics = any("Q20_pct" in d for d in seqkit_stats.values())
    has_n50 = any("N50" in d for d in seqkit_stats.values())
    has_gaps = any("sum_gap" in d for d in seqkit_stats.values())

    # General Stats Table
    general_stats_headers: Dict = {
        "num_seqs": {
            "title": "Sequences",
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
            "scale": "RdYlGn",
            "suffix": "%",
            "hidden": True,
        },
    }

    if has_quality_metrics:
        general_stats_headers["Q30_pct"] = {
            "title": "Q30%",
            "description": "Percentage of bases with quality score >= 30",
            "min": 0,
            "max": 100,
            "scale": "RdYlGn",
            "suffix": "%",
        }
        general_stats_headers["AvgQual"] = {
            "title": "Avg Qual",
            "description": "Average quality score",
            "scale": "RdYlGn",
            "hidden": True,
        }

    if has_n50:
        general_stats_headers["N50"] = {
            "title": "N50",
            "description": "N50 sequence length",
            "scale": "Oranges",
            "format": "{:,.0f}",
            "suffix": " bp",
            "hidden": True,
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
            "title": "Sequences",
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
    }

    # Add extended columns if present
    if any("Q1" in d for d in seqkit_stats.values()):
        table_headers["Q1"] = {
            "title": "Q1",
            "description": "First quartile of sequence length",
            "format": "{:,.0f}",
        }
        table_headers["Q2"] = {
            "title": "Q2 (Median)",
            "description": "Median sequence length",
            "format": "{:,.0f}",
        }
        table_headers["Q3"] = {
            "title": "Q3",
            "description": "Third quartile of sequence length",
            "format": "{:,.0f}",
        }

    if has_gaps:
        table_headers["sum_gap"] = {
            "title": "Gaps",
            "description": "Total number of gaps",
            "format": "{:,.0f}",
        }

    if has_n50:
        table_headers["N50"] = {
            "title": "N50",
            "description": "N50 sequence length",
            "format": "{:,.0f}",
        }
        if any("N50_num" in d for d in seqkit_stats.values()):
            table_headers["N50_num"] = {
                "title": "N50 Count",
                "description": "Number of sequences >= N50",
                "format": "{:,.0f}",
            }

    if has_quality_metrics:
        table_headers["Q20_pct"] = {
            "title": "Q20%",
            "description": "Percentage of bases with quality score >= 20",
            "min": 0,
            "max": 100,
            "suffix": "%",
        }
        table_headers["Q30_pct"] = {
            "title": "Q30%",
            "description": "Percentage of bases with quality score >= 30",
            "min": 0,
            "max": 100,
            "suffix": "%",
        }
        table_headers["AvgQual"] = {
            "title": "Avg Quality",
            "description": "Average quality score",
        }

    if any("GC_pct" in d for d in seqkit_stats.values()):
        table_headers["GC_pct"] = {
            "title": "GC%",
            "description": "GC content percentage",
            "min": 0,
            "max": 100,
            "suffix": "%",
        }

    if any("sum_n" in d for d in seqkit_stats.values()):
        table_headers["sum_n"] = {
            "title": "N Bases",
            "description": "Total number of N bases",
            "format": "{:,.0f}",
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


def parse_stats_report(file_content: str) -> Dict[str, Dict]:
    """
    Parse seqkit stats output file

    Returns a dictionary with sample names as keys and parsed data as values.
    Handles both tab-separated (--tabular flag) and space-separated (default) output.
    """
    parsed_data: Dict[str, Dict] = {}
    lines = file_content.strip().split("\n")

    if len(lines) < 2:
        return {}

    # Parse header line - try tab-separated first, then space-separated
    header_line = lines[0]
    if "\t" in header_line:
        headers = [h.strip() for h in header_line.split("\t")]
        use_tabs = True
    else:
        headers = header_line.split()
        use_tabs = False

    # Check if this looks like seqkit stats output
    expected_headers = ["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"]
    if not all(h in headers for h in expected_headers):
        return {}

    # Parse data lines
    for line in lines[1:]:
        if not line.strip():
            continue

        if use_tabs:
            values = line.split("\t")
        else:
            values = line.split()
        if len(values) < len(headers):
            continue

        data: Dict = {}
        sample_name = None

        for i, header in enumerate(headers):
            value = values[i].strip() if i < len(values) else ""

            if header == "file":
                # Use file name as sample name, strip path and extension
                if value == "-":
                    sample_name = "stdin"
                else:
                    # Handle both Unix and Windows path separators
                    sample_name = value.replace("\\", "/").split("/")[-1]
                # Remove common sequence file extensions
                for ext in [".fq.gz", ".fastq.gz", ".fq", ".fastq", ".fa.gz", ".fasta.gz", ".fa", ".fasta"]:
                    if sample_name.endswith(ext):
                        sample_name = sample_name[: -len(ext)]
                        break
                data["file"] = value
            elif header == "format":
                data["format"] = value
            elif header == "type":
                data["type"] = value
            elif header in [
                "num_seqs",
                "sum_len",
                "min_len",
                "max_len",
                "Q1",
                "Q2",
                "Q3",
                "sum_gap",
                "N50",
                "N50_num",
                "sum_n",
            ]:
                # Integer fields - handle comma-separated numbers
                try:
                    data[header] = int(value.replace(",", ""))
                except (ValueError, AttributeError):
                    pass
            elif header == "avg_len":
                # Float field
                try:
                    data[header] = float(value.replace(",", ""))
                except (ValueError, AttributeError):
                    pass
            elif header == "Q20(%)":
                # Quality percentage - rename for cleaner key
                try:
                    data["Q20_pct"] = float(value)
                except (ValueError, AttributeError):
                    pass
            elif header == "Q30(%)":
                try:
                    data["Q30_pct"] = float(value)
                except (ValueError, AttributeError):
                    pass
            elif header == "AvgQual":
                try:
                    data["AvgQual"] = float(value)
                except (ValueError, AttributeError):
                    pass
            elif header == "GC(%)":
                try:
                    data["GC_pct"] = float(value)
                except (ValueError, AttributeError):
                    pass

        if sample_name and data:
            parsed_data[sample_name] = data

    return parsed_data
