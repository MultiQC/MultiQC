"""MultiQC module to parse similarity matrix output by sourmash compare"""

import logging
import os
import re

import numpy

from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class CompareMixin:
    def parse_compare(self):
        """
        Modeled after vcftools relatedness2 module, which also has many samples represented in the parsed file.
        """
        matrices = {}

        for f in self.find_log_files("sourmash/compare", filehandles=True):
            labels = [x.strip() for x in f["f"]]
            if labels:
                matrix_path = re.sub(".labels.txt", "", f["f"].name)
                if not os.path.exists(matrix_path):
                    log.warning(
                        f"Found a 'labels' file expected by Sourmash: '{f['f'].name}', "
                        f"however, could not find a accompanying matrix binary file "
                        f"'{matrix_path}'. So assuming that wasn't a Sourmash result"
                    )
                    continue
                with open(matrix_path, "rb") as fh:
                    matrix = numpy.load(fh)
                # Note that "s_name" here is not a sample name, but the name of the
                # input file, that contains a comparison matrix across multiple samples.
                matrices[f["s_name"]] = (labels, matrix.tolist())
                self.add_data_source(f, section="compare")

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, f["s_name"])

        matrices = self.ignore_samples(matrices)
        if len(matrices) == 0:
            return 0

        log.info(f"Found {len(matrices)} valid compare results")

        self.write_data_file(matrices, "sourmash_compare")

        helptext = """
        Sourmash `compare` calculates the similarity score between samples. A higher score indicates a higher degree of
        similarity, up to a maximum of 1. Samples are clustered by similarity on each axis, and specific IDs can be
        found in the graph with the Highlight tab.
        """

        for name, (labels, data) in matrices.items():
            # Note that "name" here is not a sample name, but the name of the input file,
            # that contains a comparison matrix across multiple samples.
            id = name.lower().strip().replace(" ", "-").replace(".labels.txt", "")
            self.add_section(
                name=f"Sample similarity (<code>{name}</code>)",
                anchor=f"sourmash-compare-{id}",
                description=f"Heatmap of similarity values from the output of `sourmash compare` run on <code>{name}</code>",
                helptext=helptext,
                plot=heatmap.plot(
                    data,
                    xcats=labels,
                    ycats=labels,
                    pconfig={
                        "id": f"sourmash-compare-heatmap-{id}",
                        "title": "Sourmash: Compare",
                        "square": True,
                        "tt_decimals": 7,
                    },
                ),
            )

        return len(matrices)
