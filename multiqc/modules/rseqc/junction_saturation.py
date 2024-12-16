"""MultiQC submodule to parse output from RSeQC junction_saturation.py
http://rseqc.sourceforge.net/#junction-saturation-py"""

import logging
import re
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule):
    """Find RSeQC junction_saturation frequency reports and parse their data"""

    # Set up vars
    junction_saturation_all: Dict = dict()
    junction_saturation_known: Dict = dict()
    junction_saturation_novel: Dict = dict()

    # Go through files and parse data
    for f in module.find_log_files("rseqc/junction_saturation"):
        parsed = dict()
        for line in f["f"].splitlines():
            r = re.search(r"^([xyzw])=c\(([\d,]+)\)$", line)
            if r:
                parsed[r.group(1)] = [float(i) for i in r.group(2).split(",")]
        if len(parsed) == 4:
            if parsed["z"][-1] == 0:
                log.warning(f"Junction saturation data all zeroes, skipping: '{f['s_name']}'")
            else:
                if f["s_name"] in junction_saturation_all:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                module.add_data_source(f, section="junction_saturation")
                junction_saturation_all[f["s_name"]] = dict()
                junction_saturation_known[f["s_name"]] = dict()
                junction_saturation_novel[f["s_name"]] = dict()
                for k, v in enumerate(parsed["x"]):
                    junction_saturation_all[f["s_name"]][v] = parsed["z"][k]
                    junction_saturation_known[f["s_name"]][v] = parsed["y"][k]
                    junction_saturation_novel[f["s_name"]][v] = parsed["w"][k]

    # Filter to strip out ignored sample names
    junction_saturation_all = module.ignore_samples(junction_saturation_all)
    junction_saturation_known = module.ignore_samples(junction_saturation_known)
    junction_saturation_novel = module.ignore_samples(junction_saturation_novel)

    if len(junction_saturation_all) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write data to file
    module.write_data_file(junction_saturation_all, "rseqc_junction_saturation_all")
    module.write_data_file(junction_saturation_known, "junction_saturation_known")
    module.write_data_file(junction_saturation_novel, "junction_saturation_novel")

    # Add line graph to section
    pconfig = {
        "id": "rseqc_junction_saturation_plot",
        "title": "RSeQC: Junction Saturation",
        "ylab": "Number of Junctions",
        "ymin": 0,
        "xlab": "Percent of reads",
        "xmin": 0,
        "xmax": 100,
        "tt_label": "<strong>{point.x}% of reads</strong>: {point.y:.2f}",
        "data_labels": [{"name": "All Junctions"}, {"name": "Known Junctions"}, {"name": "Novel Junctions"}],
    }
    module.add_section(
        name="Junction Saturation",
        anchor="rseqc-junction_saturation",
        description="""<a href="http://rseqc.sourceforge.net/#junction-saturation-py" target="_blank">Junction Saturation</a>
            counts the number of known splicing junctions that are observed
            in each dataset. If sequencing depth is sufficient, all (annotated) splice junctions should
            be rediscovered, resulting in a curve that reaches a plateau. Missing low abundance splice
            junctions can affect downstream analysis.</p>
            <div class="alert alert-info" id="rseqc-junction_sat_single_hint">
              <span class="glyphicon glyphicon-hand-up"></span>
              Click a line to see the data side by side (as in the original RSeQC plot).
            </div><p>""",
        plot=linegraph.plot(
            [junction_saturation_all, junction_saturation_known, junction_saturation_novel],
            pconfig,
        ),
    )

    # Return number of samples found
    return len(junction_saturation_all)
