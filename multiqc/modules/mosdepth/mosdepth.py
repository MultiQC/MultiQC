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
    It also can use {prefix}.thresholds.bed.gz to generate coverage at thresholds for
    the general stats table.

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

    #chrom  start   end     region  1X      10X     20X     30X
    1       2488047 2488227 TNFRSF14        0       0       0       0
    1       2489098 2489338 TNFRSF14        0       0       0       0

    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='mosdepth', anchor='mosdepth',
            href="https://github.com/brentp/mosdepth",
            info="performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing")

        self.parse_cov_dist()

        self.all_thresholds = []
        self.parse_coverage_thresholds()
        self.general_stat_headers()

    def parse_coverage_thresholds(self):
        import pdb; pdb.set_trace()
        for f in self.find_log_files('mosdepth/thresholds'):

            print('Found ' + f)
            s_name = self.clean_s_name(f['fn'], root=None).replace('.mosdepth.thresholds.bed', '')
            hdr = list()
            parsed = list()
            for line in f['f'].split("\n"):
                if "\t" not in line:
                    continue
                if line.startswith('#'):
                    #chrom  start      end        region        0X    1X    10X   20X   30X   50X
                    hdr = line.split("\t")
                else:
                    fields = line.split("\t")
                    parsed.append(dict())
                    for h, val in zip(hdr, fields):
                        if h in ["start", "end"] or h.endswith('X'):
                            parsed[-1][h] = int(val)
            total_length = sum(d["end"] - d["start"] for d in parsed)
            thresholds = [h for h in hdr if h.endswith('X')]
            len_per_threshold = {t: sum(d[t]) for t in thresholds for d in parsed}
            for t, t_len in len_per_threshold.items():
                rate = t_len / total_length if total_length > 0 else 0
                if t not in self.all_thresholds:
                    self.all_thresholds.add(t)
                self.general_stats_data[s_name][t] = rate

            if s_name in self.general_stats_data:
                self.add_data_source(f, s_name=s_name)

    def general_stat_headers(self):
        for t in self.all_thresholds:
            self.general_stats_headers[t] = {
                'title': '&ge; {}'.format(t),
                'description': 'Fraction of genome with at least {} coverage'.format(t),
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'RdYlGn'
            }

    def parse_cov_dist(self):
        xmax = 0
        data = defaultdict(OrderedDict)
        avgdata = defaultdict(OrderedDict)

        for scope in ('region', 'global'):
            for f in self.find_log_files('mosdepth/' + scope + '_dist'):
                s_name = self.clean_s_name(f['fn'], root=None).replace('.mosdepth.' + scope + '.dist', '')
                if s_name in data:  # both region and global might exist, prioritizing region
                    continue

                for line in f['f'].split("\n"):
                    if "\t" not in line:
                        continue
                    contig, cutoff_reads, bases_fraction = line.split("\t")
                    if not contig == "total":
                        avg = avgdata[s_name].get(contig, 0) + float(bases_fraction)
                        avgdata[s_name][contig] = avg
                    y = 100.0 * float(bases_fraction)
                    x = int(cutoff_reads)
                    data[s_name][x] = y
                    if y > 1.0:
                        xmax = max(xmax, x)

                if s_name in data:
                    self.add_data_source(f, s_name=s_name)

        if data:
            self.add_section(
                name='Coverage distribution',
                anchor='mosdepth-coverage-dist',
                description='Distribution of the number of locations in the reference genome with a given depth of coverage',
                plot=linegraph.plot(data, {
                    'id': 'mosdepth-coverage-dist-id',
                    'xlab': 'Coverage (X)',
                    'ylab': '% bases in genome/regions covered by least X reads',
                    'ymax': 100,
                    'xmax': xmax,
                    'tt_label': '<b>{point.x}X</b>: {point.y}',
                })
            )
            self.add_section(
                name='Average coverage per contig',
                anchor='mosdepth-coverage-per-contig-id',
                description='Average coverage per contig or chromosome',
                plot=linegraph.plot(avgdata, {
                    'id': 'mosdepth-coverage-per-contig',
                    'xlab': 'region',
                    'ylab': 'average coverage',
                    'categories': True
                })
            )
