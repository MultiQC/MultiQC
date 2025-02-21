"""MultiQC submodule to parse output from RSeQC inner_distance.py
http://rseqc.sourceforge.net/#inner-distance-py"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC inner_distance frequency reports and parse their data"""

    inner_distance: Dict = dict()
    inner_distance_pct: Dict = dict()

    # Go through files and parse data
    for f in module.find_log_files("rseqc/inner_distance"):
        if f["s_name"] in inner_distance:
            log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
        module.add_data_source(f, section="inner_distance")
        # saving to temporary variable for SE checking later
        parsed_data = dict()
        for line in f["f"].splitlines():
            s = line.split()
            try:
                avg_pos = (float(s[0]) + float(s[1])) / 2.0
                parsed_data[avg_pos] = float(s[2])
            except Exception:
                # Don't bother running through whole file if wrong
                break

        # Only add if we actually found something i,e it was PE data
        if len(parsed_data) > 0:
            inner_distance[f["s_name"]] = parsed_data

    # Filter to strip out ignored sample names
    inner_distance = module.ignore_samples(inner_distance)

    if len(inner_distance) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None, f["s_name"])

    # Write data to file
    module.write_data_file(inner_distance, "rseqc_inner_distance")

    # Make a normalised percentage version of the data
    for s_name in inner_distance:
        inner_distance_pct[s_name] = dict()
        total = sum(inner_distance[s_name].values())
        if total == 0:
            continue
        for k, v in inner_distance[s_name].items():
            inner_distance_pct[s_name][k] = (v / total) * 100

    # Add line graph to section
    pconfig = {
        "id": "rseqc_inner_distance_plot",
        "title": "RSeQC: Inner Distance",
        "ylab": "Counts",
        "xlab": "Inner Distance (bp)",
        "tt_label": "<strong>{point.x} bp</strong>: {point.y:.2f}",
        "data_labels": [{"name": "Counts", "ylab": "Counts"}, {"name": "Percentages", "ylab": "Percentage"}],
    }
    module.add_section(
        name="Inner Distance",
        anchor="rseqc-inner_distance",
        description='<a href="http://rseqc.sourceforge.net/#inner-distance-py" target="_blank">Inner Distance</a>'
        " calculates the inner distance"
        " (or insert size) between two paired RNA reads."
        " Note that this can be negative if fragments overlap.",
        plot=linegraph.plot([inner_distance, inner_distance_pct], pconfig),
    )

    # Return number of samples found
    return len(inner_distance)
