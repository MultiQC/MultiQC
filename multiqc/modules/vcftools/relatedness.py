import csv
from collections import defaultdict
from multiqc.plots import heatmap


class Relatedness2Mixin():
    def parse_relatedness2(self):
        for f in self.find_log_files('vcftools/relatedness2', filehandles=True):
            m = _RelatednessMatrix(f)

            self.add_section(
                'Vcftools relatedness2 ' + m.name,
                'vcftools_relatedness2_' + m.name,
                plot=heatmap.plot(
                    m.matrix,
                    xcats=m.x_labels,
                    ycats=m.y_labels,
                    pconfig={'square': True, 'decimalPlaces': 3}
                ),
                description='Heatmap of RELATEDNESS_PHI values from the output of vcftools relatedness2.'
            )


class _RelatednessMatrix():
    def __init__(self, relatedness_file):
        self.name = relatedness_file['s_name']
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
