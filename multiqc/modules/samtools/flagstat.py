import logging
import re
from typing import Dict

from multiqc import config, BaseMultiqcModule
from multiqc.plots import violin

# Initialise the logger
log = logging.getLogger(__name__)


def parse_samtools_flagstat(module: BaseMultiqcModule):
    """Find Samtools flagstat logs and parse their data"""

    samtools_flagstat: Dict = dict()
    for f in module.find_log_files("samtools/flagstat"):
        parsed_data = parse_single_report(f["f"])
        if len(parsed_data) > 0:
            if f["s_name"] in samtools_flagstat:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="flagstat")
            samtools_flagstat[f["s_name"]] = parsed_data

    # Filter to strip out ignored sample names
    samtools_flagstat = module.ignore_samples(samtools_flagstat)

    if len(samtools_flagstat) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed report data to a file (restructure first)
    module.write_data_file(samtools_flagstat, "multiqc_samtools_flagstat")

    # General Stats Table
    flagstat_headers = {
        "flagstat_total": {
            "title": "Reads",
            "description": f"Total reads in the bam file ({config.read_count_desc})",
            "shared_key": "read_count",
            "hidden": True,
        },
        "mapped_passed": {
            "title": "Reads mapped",
            "description": f"Reads mapped in the bam file ({config.read_count_desc})",
            "shared_key": "read_count",
        },
        "mapped_passed_pct": {
            "title": "% Reads mapped",
            "description": "% Reads mapped in the bam file",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "hidden": True,
        },
    }
    module.general_stats_addcols(samtools_flagstat, flagstat_headers, namespace="flagstat")

    # Make a violin plot
    reads = {
        "shared_key": "read_count",
        "modify": None,
        "format": None,
        "suffix": None,
    }
    keys_counts = dict()
    keys_counts["flagstat_total"] = dict(reads, title="Total Reads")
    keys_counts["total_passed"] = dict(reads, title="Total Passed QC")
    keys_counts["mapped_passed"] = dict(reads, title="Mapped")

    if any(v.get("secondary_passed") for v in samtools_flagstat.values()):
        keys_counts["secondary_passed"] = dict(reads, title="Secondary Alignments")

    if any(v.get("supplementary_passed") for v in samtools_flagstat.values()):
        keys_counts["supplementary_passed"] = dict(reads, title="Supplementary Alignments")

    keys_counts["duplicates_passed"] = dict(reads, title="Duplicates")
    keys_counts["paired in sequencing_passed"] = dict(reads, title="Paired in Sequencing")
    keys_counts["properly paired_passed"] = dict(reads, title="Properly Paired")
    keys_counts["with itself and mate mapped_passed"] = dict(
        reads, title="Self and mate mapped", description="Reads with itself and mate mapped"
    )
    keys_counts["singletons_passed"] = dict(reads, title="Singletons")
    keys_counts["with mate mapped to a different chr_passed"] = dict(
        reads, title="Mate mapped to diff chr", description="Mate mapped to different chromosome"
    )
    keys_counts["with mate mapped to a different chr (mapQ >= 5)_passed"] = dict(
        reads, title="Diff chr (mapQ >= 5)", description="Mate mapped to different chromosome (mapQ >= 5)"
    )

    data_pct: Dict = dict()
    for sample, d in samtools_flagstat.items():
        data_pct[sample] = dict()
        total = d["flagstat_total"]
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
        name="Flagstat",
        anchor="samtools-flagstat",
        description="This module parses the output from <code>samtools flagstat</code>",
        plot=violin.plot(
            [samtools_flagstat, data_pct],
            headers=[keys_counts, keys_pct],
            pconfig={
                "id": "samtools-flagstat-dp",
                "title": "Samtools: flagstat: read count",
                "data_labels": [
                    {"name": "Read counts", "title": "Samtools flagstat: read count"},
                    {"name": "Percentage of total", "title": "Samtools flagstat: percentage of total"},
                ],
            },
        ),
    )

    # Return the number of logs that were found
    return len(samtools_flagstat)


# flagstat has one thing per line, documented here (search for flagstat):
# http://www.htslib.org/doc/samtools.html
flagstat_regexes = {
    "total": r"(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)",
    "secondary": r"(\d+) \+ (\d+) secondary",
    "supplementary": r"(\d+) \+ (\d+) supplementary",
    "duplicates": r"(\d+) \+ (\d+) duplicates",
    "mapped": r"(\d+) \+ (\d+) mapped \((.+):(.+)\)",
    "paired in sequencing": r"(\d+) \+ (\d+) paired in sequencing",
    "read1": r"(\d+) \+ (\d+) read1",
    "read2": r"(\d+) \+ (\d+) read2",
    "properly paired": r"(\d+) \+ (\d+) properly paired \((.+):(.+)\)",
    "with itself and mate mapped": r"(\d+) \+ (\d+) with itself and mate mapped",
    "singletons": r"(\d+) \+ (\d+) singletons \((.+):(.+)\)",
    "with mate mapped to a different chr": r"(\d+) \+ (\d+) with mate mapped to a different chr",
    "with mate mapped to a different chr (mapQ >= 5)": r"(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)",
}


def parse_single_report(file_obj):
    """
    Take a filename, parse the data assuming it's a flagstat file
    Returns a dictionary {'lineName_pass' : value, 'lineName_fail' : value}
    """
    parsed_data: Dict = {}

    re_groups = ["passed", "failed", "passed_pct", "failed_pct"]
    for k, r in flagstat_regexes.items():
        r_search = re.search(r, file_obj, re.MULTILINE)
        if r_search:
            for i, j in enumerate(re_groups):
                key = f"{k}_{j}"
                try:
                    val = r_search.group(i + 1).strip("% ")
                    parsed_data[key] = float(val) if ("." in val) else int(val)
                except IndexError:
                    pass  # Not all regexes have percentages
                except ValueError:
                    parsed_data[key] = float("nan")
    # Work out the total read count
    try:
        parsed_data["flagstat_total"] = parsed_data["total_passed"] + parsed_data["total_failed"]
    except KeyError:
        pass
    return parsed_data
