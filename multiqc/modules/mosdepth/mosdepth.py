""" MultiQC module to parse output from mosdepth """

from __future__ import print_function

from collections import defaultdict, OrderedDict
import logging
import fnmatch

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
from multiqc.modules.qualimap.QM_BamQC import coverage_histogram_helptext, genome_fraction_helptext
from multiqc.plots import linegraph, bargraph

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
        super(MultiqcModule, self).__init__(
            name="mosdepth",
            anchor="mosdepth",
            href="https://github.com/brentp/mosdepth",
            info="performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing",
        )

        dist_data, cov_data, xmax, perchrom_avg_data = self.parse_cov_dist()

        # Filter out any samples from --ignore-samples
        dist_data = self.ignore_samples(dist_data)
        cov_data = self.ignore_samples(cov_data)
        perchrom_avg_data = self.ignore_samples(perchrom_avg_data)

        # No samples found
        num_samples = max(len(dist_data), len(cov_data), len(perchrom_avg_data))
        if num_samples == 0:
            raise UserWarning
        log.info("Found {} reports".format(num_samples))

        if dist_data:
            self.add_section(
                name="Coverage distribution",
                anchor="mosdepth-coverage-dist",
                description="Distribution of the number of locations in the reference genome with a given depth of coverage",
                helptext=genome_fraction_helptext,
                plot=linegraph.plot(
                    dist_data,
                    {
                        "id": "mosdepth-coverage-dist-id",
                        "title": "Mosdepth: Coverage Distribution",
                        "xlab": "Coverage (X)",
                        "ylab": "% bases in genome/regions covered by at least X reads",
                        "ymax": 100,
                        "xmax": xmax,
                        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                        "smooth_points": 500,
                    },
                ),
            )
        if cov_data:
            self.add_section(
                name="Coverage plot",
                anchor="mosdepth-coverage-cov",
                description="Number of locations in the reference genome with a given depth of coverage",
                helptext=coverage_histogram_helptext,
                plot=linegraph.plot(
                    cov_data,
                    {
                        "id": "mosdepth-coverage-plot-id",
                        "title": "Mosdepth: Coverage Depth",
                        "xlab": "Coverage (X)",
                        "ylab": "% bases in genome/regions covered at X reads",
                        "ymax": 100,
                        "xmax": xmax,
                        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                        "smooth_points": 500,
                    },
                ),
            )
        if perchrom_avg_data:
            num_contigs = max([len(x.keys()) for x in perchrom_avg_data.values()])
            if num_contigs > 1:
                perchrom_plot = linegraph.plot(
                    perchrom_avg_data,
                    {
                        "id": "mosdepth-coverage-per-contig",
                        "title": "Mosdepth: Coverage per contig",
                        "xlab": "region",
                        "ylab": "average coverage",
                        "categories": True,
                        "tt_decimals": 1,
                        "tt_suffix": "x",
                        "smooth_points": 500,
                    },
                )
            else:
                perchrom_plot = bargraph.plot(
                    perchrom_avg_data,
                    pconfig={
                        "id": "mosdepth-coverage-per-contig",
                        "title": "Mosdepth: Coverage per contig",
                        "xlab": "Sample",
                        "ylab": "Average coverage",
                        "tt_suffix": "x",
                    },
                )

            self.add_section(
                name="Average coverage per contig",
                anchor="mosdepth-coverage-per-contig-id",
                description="Average coverage per contig or chromosome",
                plot=perchrom_plot,
            )
        if dist_data:
            threshs, hidden_threshs = get_cov_thresholds()
            self.genstats_cov_thresholds(dist_data, threshs, hidden_threshs)
            self.genstats_mediancov(dist_data)

    def parse_cov_dist(self):
        dist_data = defaultdict(OrderedDict)  # cumulative distribution
        cov_data = defaultdict(OrderedDict)  # absoulte (non-cumulative) coverage
        xmax = 0
        perchrom_avg_data = defaultdict(OrderedDict)  # per chromosome average coverage

        try:
            include_contigs = config.mosdepth_config.get("include_contigs", [])
            assert type(include_contigs) == list, type(include_contigs)
        except (AttributeError, TypeError, AssertionError):
            include_contigs = []
        try:
            exclude_contigs = config.mosdepth_config.get("exclude_contigs", [])
            assert type(exclude_contigs) == list, type(exclude_contigs)
        except (AttributeError, TypeError, AssertionError):
            exclude_contigs = []
        log.debug("include_contigs: {}".format(include_contigs))
        log.debug("exclude_contigs: {}".format(exclude_contigs))

        for scope in ("region", "global"):
            for f in self.find_log_files("mosdepth/" + scope + "_dist"):
                s_name = self.clean_s_name(f["fn"], f).replace(".mosdepth." + scope + ".dist", "")
                if s_name in dist_data:  # both region and global might exist, prioritizing region
                    continue

                for line in f["f"].split("\n"):
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
                        if cumcov > 1:  # require >1% to prevent long flat tail
                            xmax = max(xmax, x)
                    else:  # for per-contig plot
                        # filter out contigs based on exclusion patterns
                        if any(fnmatch.fnmatch(contig, str(pattern)) for pattern in exclude_contigs):
                            try:
                                if config.mosdepth_config["show_excluded_debug_logs"]:
                                    log.debug("Skipping excluded contig '{}'".format(contig))
                            except (AttributeError, KeyError):
                                pass
                            continue

                        # filter out contigs based on inclusion patterns
                        if len(include_contigs) > 0 and not any(
                            fnmatch.fnmatch(contig, pattern) for pattern in include_contigs
                        ):
                            # Commented out since this could be many thousands of contigs fo reach!
                            # log.debug("Skipping not included contig '{}'".format(contig))
                            continue

                        avg = perchrom_avg_data[s_name].get(contig, 0) + float(bases_fraction)
                        perchrom_avg_data[s_name][contig] = avg

                if s_name in dist_data:
                    self.add_data_source(f, s_name=s_name, section="genome_results")

        return dist_data, cov_data, xmax, perchrom_avg_data

    def genstats_cov_thresholds(self, dist_data, threshs, hidden_threshs):
        data = defaultdict(OrderedDict)
        for s_name, d in dist_data.items():
            dist_subset = {t: data for t, data in d.items() if t in threshs}
            for t in threshs:
                if int(t) in dist_subset:
                    data[s_name]["{}_x_pc".format(t)] = dist_subset[t]
                else:
                    data[s_name]["{}_x_pc".format(t)] = 0

        headers = OrderedDict()
        for t in threshs:
            headers["{}_x_pc".format(t)] = {
                "title": "&ge; {}X".format(t),
                "description": "Fraction of genome with at least {}X coverage".format(t),
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": t in hidden_threshs,
            }
        self.general_stats_addcols(data, headers)

    def genstats_mediancov(self, dist_data):
        data = defaultdict(OrderedDict)
        for s_name, d in dist_data.items():
            median_cov = None
            for this_cov, cum_pct in d.items():
                if cum_pct >= 50:
                    median_cov = this_cov
                    break
            data[s_name]["median_coverage"] = median_cov

        headers = OrderedDict()
        headers["median_coverage"] = {
            "title": "Median",
            "description": "Median coverage",
            "min": 0,
            "suffix": "X",
            "scale": "BuPu",
        }
        self.general_stats_addcols(data, headers)


def get_cov_thresholds():
    """Reads coverage thresholds from the config, otherwise sets sensible defaults"""
    try:
        threshs = config.mosdepth_config["general_stats_coverage"]
        assert type(threshs) == list
        assert len(threshs) > 0
        threshs = [int(t) for t in threshs]
        log.debug("Custom coverage thresholds: {}".format(", ".join([str(t) for t in threshs])))
    except (KeyError, AttributeError, TypeError, AssertionError):
        threshs = [1, 5, 10, 30, 50]
        log.debug("Using default coverage thresholds: {}".format(", ".join([str(t) for t in threshs])))

    try:
        hidden_threshs = config.mosdepth_config["general_stats_coverage_hidden"]
        assert type(hidden_threshs) == list
        log.debug("Hiding coverage thresholds: {}".format(", ".join([str(t) for t in hidden_threshs])))
    except (KeyError, AttributeError, TypeError, AssertionError):
        hidden_threshs = [t for t in threshs if t != 30]

    return threshs, hidden_threshs
