"""MultiQC module to parse relatedness output from vcftools relatedness"""

import csv
import logging
from collections import defaultdict

from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class Relatedness2Mixin:
    def parse_relatedness2(self):
        matrices = {}
        for f in self.find_log_files("vcftools/relatedness2", filehandles=True):
            m = _Relatedness2Matrix(f)
            if m.data and m.x_labels and m.y_labels:
                matrices[f["s_name"]] = m
            self.add_data_source(f, section="Relatedness")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        matrices = self.ignore_samples(matrices)

        if len(matrices) == 0:
            return 0

        log.info(f"Found {len(matrices)} valid relatedness2 matrices")

        # The matrices cannot be written to a file in their current format
        # self.write_data_file(matrices, "vcftools_relatedness")

        helptext = """
        `RELATEDNESS_PHI` gives a relatedness score between two samples. A higher score indicates a higher degree of
        relatedness, up to a maximum of 0.5. Samples are sorted alphabetically on each axis, and specific IDs can be
        found in the graph with the Highlight tab.
        """

        idx = 0
        for name, m in matrices.items():
            idx += 1
            self.add_section(
                name="Relatedness2",
                anchor=f"vcftools-relatedness2-{idx}",
                description="**Input:** `{}`.\n\n Heatmap of `RELATEDNESS_PHI` values from the output of vcftools relatedness2.".format(
                    name
                ),
                helptext=helptext,
                plot=heatmap.plot(
                    m.data,
                    xcats=m.x_labels,
                    ycats=m.y_labels,
                    pconfig={
                        "id": f"vcftools-relatedness2-heatmap-{idx}",
                        "title": "VCFTools: Relatedness2",
                        "square": True,
                        "tt_decimals": 7,
                    },
                ),
            )

        return len(matrices)


class _Relatedness2Matrix:
    def __init__(self, relatedness_file):
        self.data = []
        self.x_labels = set()
        self.y_labels = set()

        self.parse(relatedness_file["f"])

    def parse(self, f):
        rels = defaultdict(dict)
        r = csv.DictReader(f, delimiter="\t")
        for line in r:
            self.x_labels.add(line["INDV1"])
            self.y_labels.add(line["INDV2"])

            rels[line["INDV1"]][line["INDV2"]] = float(line["RELATEDNESS_PHI"])

        # impose alphabetical order and avoid json serialisation errors in utils.report
        self.x_labels = sorted(self.x_labels)
        self.y_labels = sorted(self.y_labels)

        for x in self.x_labels:
            line = []
            for y in self.y_labels:
                line.append(rels[x][y])
            self.data.append(line)
