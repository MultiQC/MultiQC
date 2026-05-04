"""MultiQC submodule to parse BBTools stats output (adapter/contaminant filtering)"""

import logging
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)


def parse_bbtools_stats(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools stats output files for adapter/contaminant filtering statistics.

    The stats file shows the proportion of reads that matched adapters or contaminants.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/stats", filehandles=True):
        parsed = _parse_stats_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="stats")
            data_by_sample[s_name] = parsed

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)

    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    module.add_software_version(None)

    # Create summary table with key-value metrics
    table_headers: Dict = {
        "Total": {
            "title": f"Total ({config.read_count_prefix})",
            "description": f"Total number of reads processed ({config.read_count_desc})",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "PuBu",
            "hidden": True,
        },
        "Matched": {
            "title": f"Matched ({config.read_count_prefix})",
            "description": f"Total number of reads matching adapters/contaminants ({config.read_count_desc})",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Reds",
            "hidden": True,
        },
        "Percent filtered": {
            "title": "Filtered",
            "description": "Proportion of reads filtered, matching adapters/contaminants",
            "max": 100,
            "min": 0,
            "scale": "OrRd",
            "suffix": "%",
        },
    }

    module.add_section(
        name="BBDuk Filtering Statistics",
        anchor="bbtools-stats",
        description="Proportion of reads that matched adapters/contaminants.",
        plot=table.plot(
            data_by_sample,
            table_headers,
            pconfig={
                "id": "bbtools-stats-table",
                "namespace": "BBTools",
                "title": "BBTools: Filtering Statistics",
            },
        ),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_stats")

    return len(data_by_sample)


def _parse_stats_file(f) -> Dict:
    """
    Parse a BBTools stats file.

    Expected columns: Name, Reads, ReadsPct, [Bases, BasesPct]
    Also contains key-value pairs like Total, Matched, etc.
    """
    parsed_data: Dict = {}
    adapter_data: Dict = {}

    # Expected columns in the table
    expected_cols = ["Name", "Reads", "ReadsPct"]

    for line_number, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")

        # Check if this is a header/comment line
        if parts[0].startswith("#"):
            parts[0] = parts[0][1:]  # Remove leading '#'

            # Check if it's the table header
            if parts[0] == expected_cols[0]:
                continue  # Skip header row

            # It's a key-value row
            if len(parts) == 3 and parts[0] == "File":
                # First line with three columns (paired-end) - skip
                continue
            elif len(parts) == 3:
                # Third line format: Name, Value, Percent
                try:
                    parsed_data["Percent filtered"] = float(parts[2].strip("%"))
                    parsed_data[parts[0]] = _try_parse_number(parts[1])
                except (ValueError, IndexError):
                    pass
            elif len(parts) == 2:
                # Standard key-value pair
                parsed_data[parts[0]] = _try_parse_number(parts[1])
        else:
            # Data row: Name, Reads, ReadsPct, [Bases, BasesPct]
            if len(parts) >= 3:
                name = parts[0]
                try:
                    reads = int(parts[1])
                    reads_pct = float(parts[2].strip("%"))
                    adapter_data[name] = {
                        "reads": reads,
                        "reads_pct": reads_pct,
                    }
                    if len(parts) >= 5:
                        adapter_data[name]["bases"] = int(parts[3])
                        adapter_data[name]["bases_pct"] = float(parts[4].strip("%"))
                except (ValueError, IndexError):
                    pass

    if adapter_data:
        parsed_data["adapters"] = adapter_data

    return parsed_data if parsed_data else {}


def _try_parse_number(value: str):
    """Try to parse a string as int or float, return original string if fails."""
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value
