"""MultiQC submodule to parse output from RSeQC read_duplication.py
http://rseqc.sourceforge.net/#read-duplication-py"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC read_duplication reports and parse their data"""

    read_dups: Dict = dict()

    # Go through files and parse data
    for f in module.find_log_files("rseqc/read_duplication_pos"):
        if f["f"].startswith("Occurrence	UniqReadNumber"):
            if f["s_name"] in read_dups:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="read_duplication")
            read_dups[f["s_name"]] = dict()
            for line in f["f"].splitlines():
                s = line.split()
                try:
                    if int(s[0]) <= 500:
                        read_dups[f["s_name"]][int(s[0])] = int(s[1])
                except Exception:
                    pass

    # Filter to strip out ignored sample names
    read_dups = module.ignore_samples(read_dups)

    if len(read_dups) == 0:
        return 0

    # Write data to file
    module.write_data_file(read_dups, "rseqc_read_dups")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Add line graph to section
    pconfig = {
        "smooth_points": 200,
        "id": "rseqc_read_dups_plot",
        "title": "RSeQC: Read Duplication",
        "ylab": "Number of Reads (log10)",
        "xlab": "Occurrence of read",
        "ylog": True,
        "tt_label": "<strong>{point.x} occurrences</strong>: {point.y} reads",
    }

    module.add_section(
        name="Read Duplication",
        anchor="rseqc-read_dups",
        description='<a href="http://rseqc.sourceforge.net/#read-duplication-py" target="_blank">read_duplication.py</a>'
        " calculates how many alignment positions have a certain number of exact duplicates."
        " Note - plot truncated at 500 occurrences and binned.</p>",
        plot=linegraph.plot(read_dups, pconfig),
    )

    # Return number of samples found
    return len(read_dups)
