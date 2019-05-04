""" MultiQC module to parse output from mosdepth """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
from multiqc.plots import linegraph
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Mosdepth can generate multiple outputs with a common prefix and different endings.
    The module can use first 2 (preferring "region" if exists, otherwise "global"),
    to build 2 plots: coverage distribution and per-contig average coverage.

    {prefix}.mosdepth.global.dist.txt
    a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide

    1       2       0.00
    1       1       0.00
    1       0       1.00
    total   2       0.00
    total   1       0.00
    total   0       1.00

    {prefix}.mosdepth.region.dist.txt (if --by is specified)
    same, but in regions

    1       2       0.01
    1       1       0.01
    1       0       1.00
    total   2       0.00
    total   1       0.00
    total   0       1.00

    {prefix}.per-base.bed.gz (unless -n/--no-per-base is specified)

    1       0       881481  0
    1       881481  881482  2
    1       881482  881485  4

    {prefix}.regions.bed.gz (if --by is specified)
    the mean per-region from either a BED file or windows of specified size

    1       2488047 2488227 TNFRSF14        0.00
    1       2489098 2489338 TNFRSF14        0.00

    {prefix}.quantized.bed.gz (if --quantize is specified)
    quantized output that merges adjacent bases as long as they fall in the same coverage bins e.g. (10-20)

    1       0       881481  0:1
    1       881481  881485  1:5
    1       881485  881769  5:150

    {prefix}.thresholds.bed.gz (if --thresholds is specified) - how many bases in each region are covered at the given thresholds

    #chrom  start   end     region     1X      10X     20X     30X
    1       2488047 2488227 TNFRSF14   0       0       0       0
    1       2489098 2489338 TNFRSF14   0       0       0       0

    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='mosdepth', anchor='mosdepth',
            href="https://github.com/brentp/mosdepth",
            info="performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing")

        self.parse_cov_dist()

    def parse_cov_dist(self):
        xmax = 0
        dist_data = defaultdict(OrderedDict)  # cumulative distribution
        cov_data = defaultdict(OrderedDict)  # absoulte (non-cumulative) coverage
        perchrom_avg_data = defaultdict(OrderedDict)  # per chromosome average coverage

        for scope in ('region', 'global'):
            for f in self.find_log_files('mosdepth/' + scope + '_dist'):
                s_name = self.clean_s_name(f['fn'], root=None).replace('.mosdepth.' + scope + '.dist', '')
                if s_name in dist_data:  # both region and global might exist, prioritizing region
                    continue

                for line in f['f'].split("\n"):
                    if "\t" not in line:
                        continue
                    contig, cutoff_reads, bases_fraction = line.split("\t")
                    if contig == "total":  # for global coverage distribution
                        cumcov = 100.0 * float(bases_fraction)
                        x = int(cutoff_reads)
                        dist_data[s_name][x] = cumcov
                        # converting cumulative coverage into absoulte coverage:
                        """
                        *example*              x:  cumcov:  abscov:
                        3x                     3x  0      =               0   
                        2x     -               2x  0.10   = 0.10 - 0    = 0.10  
                        1x     --------        1x  0.80   = 0.80 - 0.10 = 0.70 
                        genome ..........      0x  1.00   = 1.00 - 0.80 = 0.20
                        """
                        if x + 1 not in dist_data[s_name]:
                            cov_data[s_name][x] = cumcov
                        else:
                            cov_data[s_name][x] = cumcov - dist_data[s_name][x + 1]
                        if cumcov > 1.0:
                            xmax = max(xmax, x)
                    else:  # for per-contig plot
                        avg = perchrom_avg_data[s_name].get(contig, 0) + float(bases_fraction)
                        perchrom_avg_data[s_name][contig] = avg

                if s_name in dist_data:
                    self.add_data_source(f, s_name=s_name)

        if dist_data:
            self.add_section(
                name='Coverage distribution',
                anchor='mosdepth-coverage-dist',
                description='Distribution of the number of locations in the reference genome with a given depth of coverage',
                plot=linegraph.plot(dist_data, {
                    'id': 'mosdepth-coverage-dist-id',
                    'xlab': 'Coverage (X)',
                    'ylab': '% bases in genome/regions covered by least X reads',
                    'ymax': 100,
                    'xmax': xmax,
                    'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
                })
            )
            self.add_section(
                name='Coverage plot',
                anchor='mosdepth-coverage-cov',
                description='Number of locations in the reference genome with a given depth of coverage',
                plot=linegraph.plot(cov_data, {
                    'id': 'mosdepth-coverage-plot-id',
                    'xlab': 'Coverage (X)',
                    'ylab': '% bases in genome/regions covered at X reads',
                    'ymax': 100,
                    'xmax': xmax,
                    'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
                })
            )
            self.add_section(
                name='Average coverage per contig',
                anchor='mosdepth-coverage-per-contig-id',
                description='Average coverage per contig or chromosome',
                plot=linegraph.plot(perchrom_avg_data, {
                    'id': 'mosdepth-coverage-per-contig',
                    'xlab': 'region',
                    'ylab': 'average coverage',
                    'categories': True,
                    'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
                })
            )
