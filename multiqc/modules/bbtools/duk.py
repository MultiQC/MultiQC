"""MultiQC submodule to parse output from BBDuk"""

import logging
import re
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"Version ([\d\.]+)"


def parse_bbtools_duk(module: BaseMultiqcModule) -> int:
    """
    Parse BBDuk stdout logs and extract filtering statistics.

    BBDuk (Decontamination Using Kmers) combines common data-quality-related
    trimming, filtering, and masking operations into a single high-performance tool.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/duk", filehandles=True):
        parsed_data = _parse_duk_log(f["f"], f["s_name"], module, f)
        if parsed_data:
            s_name = parsed_data.pop("_s_name")
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="duk")
            data_by_sample[s_name] = parsed_data

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)

    if len(data_by_sample) == 0:
        return 0

    # General Stats Table
    general_stats_headers: Dict = {
        "Total Removed bases percent": {
            "title": "Bases Removed (%)",
            "description": "Percentage of bases removed after filtering",
            "scale": "YlOrBr",
            "max": 100,
            "suffix": "%",
        },
        "Total Removed bases": {
            "title": f"Bases Removed ({config.base_count_prefix})",
            "description": f"Total bases removed ({config.base_count_desc})",
            "scale": "Reds",
            "shared_key": "base_count",
            "modify": lambda x: x * config.base_count_multiplier,
            "hidden": True,
        },
        "Total Removed reads percent": {
            "title": "Reads Removed (%)",
            "description": "Percentage of reads removed after filtering",
            "scale": "OrRd",
            "max": 100,
            "suffix": "%",
        },
        "Total Removed reads": {
            "title": f"Reads Removed ({config.read_count_prefix})",
            "description": f"Total reads removed ({config.read_count_desc})",
            "scale": "Reds",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        },
        "Input reads": {
            "title": f"Input Reads ({config.read_count_prefix})",
            "description": f"Total number of input reads to BBDuk ({config.read_count_desc})",
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "hidden": True,
        },
    }

    general_stats_headers = module.get_general_stats_headers(all_headers=general_stats_headers)
    if general_stats_headers:
        module.general_stats_addcols(data_by_sample, general_stats_headers, namespace="duk")

    # Bar graph of filtered reads
    cats = [
        "Result",
        "QTrimmed",
        "KTrimmed",
        "Trimmed by overlap",
        "Low quality discards",
        "Low entropy discards",
    ]

    pconfig = {
        "id": "bbtools-duk-filtered-barplot",
        "title": "BBTools: BBDuk Filtered Reads",
        "ylab": "Number of Reads",
        "data_labels": [
            {"name": "Reads", "ylab": "Number of Reads"},
            {"name": "Bases", "ylab": "Number of Base Pairs"},
        ],
    }

    module.add_section(
        name="BBDuk: Filtered Reads",
        anchor="bbtools-duk",
        description="The number of reads removed by various BBDuk filters.",
        plot=bargraph.plot(
            [data_by_sample, data_by_sample],
            [
                [f"{cat} reads" for cat in cats],
                [f"{cat} bases" for cat in cats],
            ],
            pconfig,
        ),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_duk")

    return len(data_by_sample)


def _parse_duk_log(f, s_name: str, module: BaseMultiqcModule, file_info) -> Dict:
    """Parse a BBDuk stdout log file."""
    parsed_data: Dict = {}
    current_s_name = s_name

    for line in f:
        # Extract sample name from command line if available
        if "jgi.BBDuk" in line and "in1=" in line:
            current_s_name = line.split("in1=")[1].split(" ")[0]
            current_s_name = module.clean_s_name(current_s_name, file_info)

        # Extract version
        if line.startswith("Version"):
            version_match = re.search(VERSION_REGEX, line)
            if version_match:
                module.add_software_version(version_match.group(1), current_s_name)

        # Parse input reads/bases
        if "Input:" in line:
            matches = re.search(r"Input:\s+(\d+) reads\s+(\d+) bases", line)
            if matches:
                parsed_data["Input reads"] = int(matches.group(1))
                parsed_data["Input bases"] = int(matches.group(2))

        # Parse filtering categories (only after we've seen Input)
        elif "Input reads" in parsed_data:
            cats = [
                "QTrimmed",
                "KTrimmed",
                "Trimmed by overlap",
                "Low quality discards",
                "Low entropy discards",
                "Total Removed",
                "Result",
            ]
            for cat in cats:
                matches = re.search(rf"{cat}:\s+(\d+) reads \(([\d\.]+)%\)\s+(\d+) bases \(([\d\.]+)%\)", line)
                if matches:
                    parsed_data[f"{cat} reads"] = int(matches.group(1))
                    parsed_data[f"{cat} reads percent"] = float(matches.group(2))
                    parsed_data[f"{cat} bases"] = int(matches.group(3))
                    parsed_data[f"{cat} bases percent"] = float(matches.group(4))
                    break

        # Stop parsing at end of BBDuk output
        elif "Reads Processed:" in line:
            break

    if parsed_data:
        parsed_data["_s_name"] = current_s_name

    return parsed_data
