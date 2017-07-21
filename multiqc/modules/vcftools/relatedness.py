#!/usr/bin/env python

""" MultiQC module to parse relatedness output from vcftools relatedness """

import csv
import logging
from collections import defaultdict
from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class Relatedness2Mixin():
    def parse_relatedness2(self):
        matrices = {}
        for f in self.find_log_files('vcftools/relatedness2', filehandles=True):
            m = _RelatednessMatrix(f)
            if m.matrix and m.x_labels and m.y_labels:
                matrices[f['s_name']] = m

        matrices = self.ignore_samples(matrices)
        log.info('Found %s relatedness2 files', len(matrices))

        helptext = '''
        RELATEDNESS_PHI gives a relatedness score between two samples. A higher score indicates a higher degree of
        relatedness, up to a maximum of 0.5.
        '''

        for name, m in matrices.items():
            self.add_section(
                name='Vcftools relatedness2 ' + name,
                anchor='vcftools_relatedness2_' + name,
                description='Heatmap of RELATEDNESS_PHI values from the output of vcftools relatedness2.',
                helptext=helptext,
                plot=heatmap.plot(
                    m.matrix,
                    xcats=m.x_labels,
                    ycats=m.y_labels,
                    pconfig={'square': True, 'decimalPlaces': 3}
                )
            )

        return len(matrices)


class _RelatednessMatrix():
    def __init__(self, relatedness_file):
        self.matrix = []
        self.x_labels = set()
        self.y_labels = set()

        self.parse(relatedness_file['f'])

    def parse(self, f):
        rels = defaultdict(dict)
        r = csv.DictReader(f, delimiter='\t')
        for line in r:
            self.x_labels.add(line['INDV1'])
            self.y_labels.add(line['INDV2'])

            rels[line['INDV1']][line['INDV2']] = float(line['RELATEDNESS_PHI'])

        self.x_labels = sorted(self.x_labels)
        self.y_labels = sorted(self.y_labels)
        for x in self.x_labels:
            line = []
            for y in self.y_labels:
                line.append(rels[x][y])
            self.matrix.append(line)
