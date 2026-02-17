import logging
import re
from typing import Any, Dict, Optional

from multiqc import BaseMultiqcModule, config
from multiqc.plots import violin

# Initialise the logger
log = logging.getLogger(__name__)

# ampliconclip has value thing per line, documented here (search for ampliconclip):
# http://www.htslib.org/doc/samtools.html
ampliconclip_headers = {
    "ampliconclip_total_reads": {
        "title": "Total Reads",
        "description": f"Total Reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "TOTAL READS",
    },
    "ampliconclip_total_clipped": {
        "title": "Total Clipped",
        "description": "Total Clipped ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "TOTAL CLIPPED",
    },
    "ampliconclip_forward_clipped": {
        "title": "Forward Clipped",
        "description": "Forward Clipped ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "FORWARD CLIPPED",
    },
    "ampliconclip_reverse_clipped": {
        "title": "Reverse Clipped",
        "description": "Reverse Clipped ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "REVERSE CLIPPED",
    },
    "ampliconclip_both_clipped": {
        "title": "Both Clipped",
        "description": "Both Clipped ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "BOTH CLIPPED",
    },
    "ampliconclip_not_clipped": {
        "title": "Not Clipped",
        "description": "Not Clipped ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "NOT CLIPPED",
    },
    "ampliconclip_excluded_reads": {
        "title": "Excluded Reads",
        "description": "Excluded Reads({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "EXCLUDED",
    },
    "ampliconclip_filtered_reads": {
        "title": "Filtered Reads",
        "description": "Filtered Reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "FILTERED",
    },
    "ampliconclip_failed_reads": {
        "title": "Failed Reads",
        "description": "Failed Reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "FAILED",
    },
    "ampliconclip_written_reads": {
        "title": "Written Reads",
        "description": "Written Reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "hidden": True,
        "source_col": "WRITTEN",
    },
}


def parse_samtools_ampliconclip(module: BaseMultiqcModule):
    """Find Samtools ampliconclip logs and parse their data"""

    samtools_ampliconclip: Dict = dict()
    for f in module.find_log_files("samtools/ampliconclip"):
        parsed_data = parse_single_report(f["f"])
        if len(parsed_data) > 0:
            if f["s_name"] in samtools_ampliconclip:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="ampliconclip")
            samtools_ampliconclip[f["s_name"]] = parsed_data

    # Filter to strip out ignored sample names
    samtools_ampliconclip = module.ignore_samples(samtools_ampliconclip)

    if len(samtools_ampliconclip) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # General Stats Table

    # Get general stats headers using the utility function, will read config.general_stats_columns
    general_stats_headers = module.get_general_stats_headers(
        all_headers={
            data_dict["source_col"]: {
                "title": data_dict["title"],
                "description": data_dict["description"],
                "shared_key": data_dict["shared_key"],
                "hidden": data_dict["hidden"],
            }
            for data_key, data_dict in ampliconclip_headers.items()
        }
    )

    # Add headers to general stats table
    if general_stats_headers:
        module.general_stats_addcols(samtools_ampliconclip, general_stats_headers, namespace="ampliconclip")

    # Make a violin plot
    keys_counts: Dict[str, Dict[str, Any]] = {
        key: {
            "shared_key": "read_count",
            "modify": None,
            "format": None,
            "suffix": None,
            "title": data["title"],
        }
        for key, data in ampliconclip_headers.items()
    }

    data_pct: Dict = dict()
    for sample, d in samtools_ampliconclip.items():
        data_pct[sample] = dict()
        total = d["ampliconclip_total_reads"]
        if total > 0:
            for metric, cnt in d.items():
                data_pct[sample][f"{metric}_pct"] = cnt / total * 100
    keys_pct = {
        f"{metric}_pct": dict(
            title=header["title"],
            min=0,
            max=100,
            suffix="%",
            shared_key=None,
        )
        for metric, header in keys_counts.items()
    }

    module.add_section(
        name="Amplicon Clip",
        anchor="samtools-ampliconclip",
        description="This module parses the output from <code>samtools ampliconclip</code>",
        plot=violin.plot(
            samtools_ampliconclip,
            headers=keys_counts,
            pconfig={
                "id": "samtools-ampliconclip-table",
                "title": "Samtools: ampliconclip: read count",
            },
        ),
    )

    module.add_section(
        name="Amplicon Clip: Percentage of total",
        anchor="samtools-ampliconclip-pct",
        description="This module parses the output from <code>samtools ampliconclip</code>",
        plot=violin.plot(
            data_pct,
            headers=keys_pct,
            pconfig={
                "id": "samtools-ampliconclip-pct-table",
                "title": "Samtools: ampliconclip: percentage of total",
            },
        ),
    )

    # Write parsed report data to a file (restructure first)
    module.write_data_file(samtools_ampliconclip, "multiqc_samtools_ampliconclip")

    # Return the number of logs that were found
    return len(samtools_ampliconclip)


def parse_single_report(file_obj):
    """
    Take a filename, parse the data assuming it's a ampliconclip file
    Returns a dictionary of metrics to value
    """
    parsed_data: Dict[str, Optional[int]] = {k: None for k in ampliconclip_headers}

    source_to_key = {data_dict["source_col"]: data_key for data_key, data_dict in ampliconclip_headers.items()}
    source_to_key["COMMAND"] = None

    for line in file_obj.splitlines():
        source, value = line.split(":", maxsplit=1)
        data_key = source_to_key[source]
        if data_key is None:
            continue
        parsed_data[data_key] = int(value)

    # check that all the values were found
    for k, v in parsed_data.items():
        if v is None:
            raise ValueError(f"Missing value '{k}' for samtools ampliconclip")

    return parsed_data
