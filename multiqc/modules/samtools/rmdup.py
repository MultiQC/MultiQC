import logging
import re
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_samtools_rmdup(module: BaseMultiqcModule):
    """Find Samtools rmdup logs and parse their data"""

    samtools_rmdup: Dict = dict()
    for f in module.find_log_files("samtools/rmdup", filehandles=True):
        # Example below:
        # [bam_rmdupse_core] 26602816 / 103563641 = 0.2569 in library '   '
        dups_regex = r"\[bam_rmdups?e?_core\] (\d+) / (\d+) = (\d+\.\d+) in library '(.*)'"
        s_name = f["s_name"]
        for line in f["f"]:
            match = re.search(dups_regex, line)
            if match:
                library_name = match.group(4).strip()
                if library_name != "":
                    s_name = library_name
                if s_name in samtools_rmdup:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                module.add_data_source(f, s_name)
                samtools_rmdup[s_name] = {}
                samtools_rmdup[s_name]["n_dups"] = int(match.group(1))
                samtools_rmdup[s_name]["n_tot"] = int(match.group(2))
                samtools_rmdup[s_name]["n_unique"] = int(match.group(2)) - int(match.group(1))
                samtools_rmdup[s_name]["pct_dups"] = float(match.group(3)) * 100

    # Filter to strip out ignored sample names
    samtools_rmdup = module.ignore_samples(samtools_rmdup)

    if len(samtools_rmdup) == 0:
        return 0

    # Write parsed report data to a file
    module.write_data_file(samtools_rmdup, "multiqc_samtools_rmdup")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Make a bar plot showing duplicates
    keys = {
        "n_unique": {"name": "Non-duplicated reads"},
        "n_dups": {"name": "Duplicated reads"},
    }
    pconfig = {
        "id": "samtools_rmdup_plot",
        "title": "Samtools: rmdup: Duplicate alignments",
        "ylab": "Number of reads",
        "y_decimals": False,
    }
    module.add_section(
        name="Duplicates removed",
        anchor="samtools-rmdup",
        plot=bargraph.plot(samtools_rmdup, keys, pconfig),
    )

    # Add a column to the General Stats table
    # General Stats Table
    stats_headers = {
        "pct_dups": {
            "title": "Duplicates",
            "description": "Percent of duplicate alignments",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "OrRd",
        }
    }
    module.general_stats_addcols(samtools_rmdup, stats_headers, namespace="rmdup")

    return len(samtools_rmdup)
