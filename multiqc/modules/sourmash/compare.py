#!/usr/bin/env python

""" MultiQC module to parse similarity matrix output by sourmash compare """

import logging
import numpy
import re

from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)

class compare:
    def parse_compare(self):
        """
        Modeled after vcftools relatedness2 module, which also has many samples represented.
        The alternative would be to have the root of the labels.txt file be ignored
        """
        matrices = {}
        
        for f in self.find_log_files("sourmash/compare", filehandles=True):
            m = compare2matrix(f)
            if m.data and m.comparelabels:
                 matrices[f["s_name"]] = m
            self.add_data_source(f, section="compare")

        matrices = self.ignore_samples(matrices)

        if len(matrices) == 0:
            return 0

        log.info("Found {} valid compare matrices".format(len(matrices)))

        # The matrices cannot be written to a file in their current format
        # self.write_data_file(matrices, "sourmash_compare")

        helptext = """
        Sourmash compare outputs a similarity score between two samples. A higher score indicates a higher degree of
        similarity, up to a maximum of 1. Samples are clustered by similarity on each axis, and specific IDs can be
        found in the graph with the Highlight tab.
        """

        idx = 0
        for name, m in matrices.items():
            idx += 1
            self.add_section(
                name = "compare",
                anchor = "sourmash-compare-{}".format(idx),
                description = "**Input:** `{}`.\n\n Heatmap of similarity values from the output of sourmash compare".format(
                    name
                ),
                helptext = helptext,
                plot = heatmap.plot(
                     m.data,
                     xcats = m.comparelabels,
                     ycats = m.comparelabels,
                     pconfig = {
                        "id": "sourmash-compare-heatmap-{}".format(idx),
                        "title": "sourmash: compare",
                        "square": True,
                        "decimalPlaces": 7,
                     },
                ),
            )

        return len(matrices)


class compare2matrix:
    def __init__(self, compare_file):
        self.data = []
        self.comparelabels = set()

        self.load_matrix_and_labels(compare_file["f"])

    def load_matrix_and_labels(self, f):
        """
        source for first two lines: https://github.com/sourmash-bio/sourmash/blob/9083d20aabcb77c67ba050b727efdd3f5d0a0398/src/sourmash/fig.py#L13
        """
        self.comparelabels = [x.strip() for x in f] 
        basefile = re.sub(".labels.txt", "", str(f.name))
        comparematrix = numpy.load(open(basefile, 'rb'))
        comparedict = {}
        for i in range(len(self.comparelabels)):
            comparevalues = list(comparematrix[i])
            res = {self.comparelabels[i]: comparevalues[i] for i in range(len(self.comparelabels))}
            comparedict[self.comparelabels[i]] = res
    
        # impose alphabetical order and avoid json serialisation errors in utils.report
        self.comparelabels = sorted(self.comparelabels)

        for x in self.comparelabels:
            line = []
            for y in self.comparelabels:
                line.append(comparedict[x][y])
            self.data.append(line)
