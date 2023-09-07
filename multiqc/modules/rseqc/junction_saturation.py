""" MultiQC submodule to parse output from RSeQC junction_saturation.py
http://rseqc.sourceforge.net/#junction-saturation-py """

import logging
import re
from collections import OrderedDict

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find RSeQC junction_saturation frequency reports and parse their data"""

    # Set up vars
    self.junction_saturation_all = dict()
    self.junction_saturation_known = dict()
    self.junction_saturation_novel = dict()

    # Go through files and parse data
    for f in self.find_log_files("rseqc/junction_saturation"):
        parsed = dict()
        for l in f["f"].splitlines():
            r = re.search(r"^([xyzw])=c\(([\d,]+)\)$", l)
            if r:
                parsed[r.group(1)] = [float(i) for i in r.group(2).split(",")]
        if len(parsed) == 4:
            if parsed["z"][-1] == 0:
                log.warning("Junction saturation data all zeroes, skipping: '{}'".format(f["s_name"]))
            else:
                if f["s_name"] in self.junction_saturation_all:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f, section="junction_saturation")
                self.junction_saturation_all[f["s_name"]] = OrderedDict()
                self.junction_saturation_known[f["s_name"]] = OrderedDict()
                self.junction_saturation_novel[f["s_name"]] = OrderedDict()
                for k, v in enumerate(parsed["x"]):
                    self.junction_saturation_all[f["s_name"]][v] = parsed["z"][k]
                    self.junction_saturation_known[f["s_name"]][v] = parsed["y"][k]
                    self.junction_saturation_novel[f["s_name"]][v] = parsed["w"][k]

    # Filter to strip out ignored sample names
    self.junction_saturation_all = self.ignore_samples(self.junction_saturation_all)
    self.junction_saturation_known = self.ignore_samples(self.junction_saturation_known)
    self.junction_saturation_novel = self.ignore_samples(self.junction_saturation_novel)

    if len(self.junction_saturation_all) > 0:
        # Write data to file
        self.write_data_file(self.junction_saturation_all, "rseqc_junction_saturation_all")
        self.write_data_file(self.junction_saturation_known, "junction_saturation_known")
        self.write_data_file(self.junction_saturation_novel, "junction_saturation_novel")

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
            "cursor": "pointer",
            "click_func": "single_sample_plot",
        }
        self.add_section(
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
                [self.junction_saturation_all, self.junction_saturation_known, self.junction_saturation_novel], pconfig
            ),
        )

    # Return number of samples found
    return len(self.junction_saturation_all)
