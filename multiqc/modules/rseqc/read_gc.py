"""MultiQC submodule to parse output from RSeQC read_GC.py
http://rseqc.sourceforge.net/#read-gc-py"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC read_GC reports and parse their data"""

    # Set up vars
    read_gc: Dict = dict()
    read_gc_pct: Dict = dict()

    # Go through files and parse data
    for f in module.find_log_files("rseqc/read_gc"):
        if f["f"].startswith("GC%	read_count"):
            gc = list()
            counts = list()
            for line in f["f"].splitlines():
                s = line.split()
                try:
                    gc.append(float(s[0]))
                    counts.append(float(s[1]))
                except Exception:
                    pass
            if len(gc) > 0:
                sorted_gc_keys = sorted(range(len(gc)), key=lambda k: gc[k])
                total = sum(counts)
                if f["s_name"] in read_gc:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                module.add_data_source(f, section="read_GC")
                read_gc[f["s_name"]] = dict()
                read_gc_pct[f["s_name"]] = dict()
                for i in sorted_gc_keys:
                    read_gc[f["s_name"]][gc[i]] = counts[i]
                    read_gc_pct[f["s_name"]][gc[i]] = (counts[i] / total) * 100

    # Filter to strip out ignored sample names
    read_gc = module.ignore_samples(read_gc)

    if len(read_gc) == 0:
        return 0

    # Write data to file
    module.write_data_file(read_gc, "rseqc_read_gc")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Add line graph to section
    pconfig = {
        "id": "rseqc_read_gc_plot",
        "title": "RSeQC: Read GC Content",
        "ylab": "Number of Reads",
        "xlab": "GC content (%)",
        "xmin": 0,
        "xmax": 100,
        "tt_label": "<strong>{point.x}% GC</strong>: {point.y:.2f}",
        "data_labels": [
            {"name": "Counts", "ylab": "Number of Reads"},
            {"name": "Percentages", "ylab": "Percentage of Reads"},
        ],
    }
    module.add_section(
        name="Read GC Content",
        anchor="rseqc-read_gc",
        description='<a href="http://rseqc.sourceforge.net/#read-gc-py" target="_blank">read_GC</a>'
        " calculates a histogram of read GC content.</p>",
        plot=linegraph.plot([read_gc, read_gc_pct], pconfig),
    )

    # Return number of samples found
    return len(read_gc)
