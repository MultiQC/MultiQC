#!/usr/bin/env python

""" MultiQC module to parse similarity matrix output by sourmash compare """

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
                    raise RuntimeError(
                        f"Expected to find the matrix binary file {matrix_path} "
                        f"complementing the labels file {f['f'].name}"
                    )
                with open(matrix_path, "rb") as fh:
                    matrix = numpy.load(fh)
                matrices[f["s_name"]] = (labels, matrix.tolist())
                self.add_data_source(f, section="compare")

        matrices = self.ignore_samples(matrices)
        if len(matrices) == 0:
            return 0

        log.info("Found {} valid compare results".format(len(matrices)))

        self.write_data_file(matrices, "sourmash_compare")

        helptext = """
        Sourmash compare outputs a similarity score between two samples. A higher score indicates a higher degree of
        similarity, up to a maximum of 1. Samples are clustered by similarity on each axis, and specific IDs can be
        found in the graph with the Highlight tab.
        """

        idx = 0
        for name, (labels, data) in matrices.items():
            idx += 1
            self.add_section(
                name="Compare: Sample Similarity",
                anchor="sourmash-compare-{}".format(idx),
                description=f"**Input:** `{name}`.\n\n Heatmap of similarity values from the output of sourmash compare",
                helptext=helptext,
                plot=heatmap.plot(
                    data,
                    xcats=labels,
                    ycats=labels,
                    pconfig={
                        "id": "sourmash-compare-heatmap-{}".format(idx),
                        "title": "Sourmash: Compare",
                        "square": True,
                        "decimalPlaces": 7,
                    },
                ),
            )

        return len(matrices)
