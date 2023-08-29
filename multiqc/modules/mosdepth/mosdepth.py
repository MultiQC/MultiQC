""" MultiQC module to parse output from mosdepth """


import fnmatch
import logging
from collections import OrderedDict, defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
from multiqc.modules.qualimap.QM_BamQC import coverage_histogram_helptext, genome_fraction_helptext
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)


def read_config():
    cfg = getattr(config, "mosdepth_config", dict())
    if type(cfg) != dict:
        return {}

    cfg["include_contigs"] = cfg.get("include_contigs", [])
    if type(cfg["include_contigs"]) != list:
        cfg["include_contigs"] = []

    cfg["exclude_contigs"] = cfg.get("exclude_contigs", [])
    if type(cfg["exclude_contigs"]) != list:
        cfg["exclude_contigs"] = []

    cfg["xchr"] = cfg.get("xchr", None)
    if type(cfg["xchr"]) != str:
        cfg["xchr"] = None

    cfg["ychr"] = cfg.get("ychr", None)
    if type(cfg["ychr"]) != str:
        cfg["ychr"] = None

    if cfg["include_contigs"]:
        log.debug("Trying to include these contigs in mosdepth: {}".format(", ".join(cfg["include_contigs"])))
    if cfg["exclude_contigs"]:
        log.debug("Excluding these contigs from mosdepth: {}".format(", ".join(cfg["exclude_contigs"])))
    if cfg["xchr"]:
        log.debug('Using "{}" as X chromosome name'.format(cfg["xchr"]))
    if cfg["ychr"]:
        log.debug('Using "{}" as Y chromosome name'.format(cfg["ychr"]))

    cutoff = cfg.get("perchrom_fraction_cutoff", 0.0)
    try:
        cutoff = float(cutoff)
    except ValueError:
        cutoff = 0.0
    if cutoff != 0.0:
        log.debug(f"Setting mosdepth coverage cutoff to display the contigs to " f"{cutoff * 100.0}%")
    cfg["perchrom_fraction_cutoff"] = cutoff

    return cfg


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
            doi="10.1093/bioinformatics/btx699",
        )

        self.cfg = read_config()
        genstats_headers = defaultdict(OrderedDict)
        genstats = defaultdict(OrderedDict)  # mean coverage

        # Parse mean coverage
        for f in self.find_log_files("mosdepth/summary"):
            s_name = self.clean_s_name(f["fn"], f)
            for line in f["f"].splitlines():
                if line.startswith("total\t"):
                    contig, length, bases, mean, min_cov, max_cov = line.split("\t")
                    genstats[s_name]["mean_coverage"] = mean
                    self.add_data_source(f, s_name=s_name, section="summary")
        # Filter out any samples from --ignore-samples
        genstats = defaultdict(OrderedDict, self.ignore_samples(genstats))
        samples_found = set(genstats.keys())

        data_by_scope = {}
        for scope in "global", "region":
            data_by_scope[scope] = self.parse_cov_dist(scope)
            data_by_scope[scope] = [defaultdict(OrderedDict, self.ignore_samples(d)) for d in data_by_scope[scope]]
            for d in data_by_scope[scope]:
                samples_found.update(set(d.keys()))

        # No samples found
        if len(samples_found) == 0:
            raise UserWarning
        log.info(f"Found {len(samples_found)} reports")

        for scope, data in data_by_scope.items():
            if scope == "global":
                descr_suf = ". Calculated across the entire genome length"
                title_suf = ""
                id_suf = ""
                fn_suf = ""
            else:
                descr_suf = ". Calculated across the target regions"
                title_suf = " (regions only)"
                id_suf = "-regions"
                fn_suf = "_regions"

            cumcov_dist_data, cov_dist_data, perchrom_avg_data, xy_cov = data
            if cumcov_dist_data:
                xmax = 0
                for sample, data in cumcov_dist_data.items():
                    for x, cumcov in data.items():
                        if cumcov > 1:  # require >1% to prevent long flat tail
                            xmax = max(xmax, x)

                # Write data to file, sort columns numerically and convert to strings
                cumcov_dist_data_writeable = {
                    sample: {str(k): v for k, v in sorted(data.items())} for sample, data in cumcov_dist_data.items()
                }
                self.write_data_file(cumcov_dist_data_writeable, f"mosdepth_cumcov{fn_suf}_dist")

                self.add_section(
                    name=f"Cumulative coverage distribution{title_suf}",
                    anchor=f"mosdepth-cumcoverage{id_suf}-dist",
                    description=(
                        f"Proportion of bases in the reference genome with, "
                        f"at least, a given depth of coverage{descr_suf}"
                    ),
                    helptext=genome_fraction_helptext,
                    plot=linegraph.plot(
                        cumcov_dist_data,
                        {
                            "id": f"mosdepth-cumcoverage-dist{id_suf}-id",
                            "title": f"Mosdepth: Cumulative coverage distribution{title_suf}",
                            "xlab": "Cumulative Coverage (X)",
                            "ylab": "% bases in genome/regions covered by at least X reads",
                            "ymax": 100,
                            "xmax": xmax,
                            "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                            "smooth_points": 500,
                        },
                    ),
                )

                assert cov_dist_data, "cov_dist_data is built from the same source and must exist here"
                # Write data to file, sort columns numerically and convert to strings
                cov_dist_data_writeable = {
                    sample: {str(k): v for k, v in sorted(data.items())} for sample, data in cov_dist_data.items()
                }
                self.write_data_file(cov_dist_data_writeable, f"mosdepth_cov{fn_suf}_dist")

                # Set ymax so that zero coverage values are ignored.
                ymax = 0
                for data in cov_dist_data.values():
                    positive_cov = [percent for cov, percent in data.items() if cov > 0]
                    if positive_cov:
                        ymax = max(ymax, max(positive_cov))

                self.add_section(
                    name=f"Coverage distribution{title_suf}",
                    anchor=f"mosdepth-coverage-dist{id_suf}-cov",
                    description=(
                        f"Proportion of bases in the reference genome with a given " f"depth of coverage{descr_suf}"
                    ),
                    helptext=coverage_histogram_helptext,
                    plot=linegraph.plot(
                        cov_dist_data,
                        {
                            "id": f"mosdepth-coverage-dist{id_suf}-id",
                            "title": f"Mosdepth: Coverage distribution{title_suf}",
                            "xlab": "Coverage (X)",
                            "ylab": "% bases in genome/regions covered by X reads",
                            "ymax": ymax * 1.05,
                            "yCeiling": 100,
                            "xmax": xmax,
                            "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                            "smooth_points": 500,
                        },
                    ),
                )
            if perchrom_avg_data:
                # Write data to file
                self.write_data_file(perchrom_avg_data, f"mosdepth_perchrom{fn_suf}")

                num_contigs = max([len(x.keys()) for x in perchrom_avg_data.values()])
                if num_contigs > 1:
                    perchrom_plot = linegraph.plot(
                        perchrom_avg_data,
                        {
                            "id": f"mosdepth-coverage-per-contig{id_suf}",
                            "title": f"Mosdepth: Coverage per contig{title_suf}",
                            "xlab": "Region",
                            "ylab": "Average Coverage",
                            "categories": True,
                            "tt_decimals": 1,
                            "tt_suffix": "x",
                            "smooth_points": 500,
                            "logswitch": True,
                            "hide_zero_cats": False,
                        },
                    )
                else:
                    perchrom_plot = bargraph.plot(
                        perchrom_avg_data,
                        pconfig={
                            "id": f"mosdepth-coverage-per-contig{id_suf}",
                            "title": f"Mosdepth: Coverage per contig{title_suf}",
                            "xlab": "Sample",
                            "ylab": "Average Coverage",
                            "tt_suffix": "x",
                            "hide_zero_cats": False,
                        },
                    )

                self.add_section(
                    name=f"Average coverage per contig{title_suf}",
                    anchor=f"mosdepth-coverage-per-contig{id_suf}-id",
                    description=f"Average coverage per contig or chromosome{id_suf}",
                    plot=perchrom_plot,
                )

            if xy_cov:
                xy_keys = OrderedDict()
                xy_keys["x"] = {"name": self.cfg.get("xchr", "Chromosome X")}
                xy_keys["y"] = {"name": self.cfg.get("xchr", "Chromosome Y")}
                pconfig = {
                    "id": f"mosdepth-xy-coverage-plot{id_suf}",
                    "title": f"Mosdepth: chrXY coverage{title_suf}",
                    "ylab": f"Percent of X+Y coverage",
                    "cpswitch_counts_label": "Coverage",
                    "cpswitch_percent_label": "Percent of X+Y coverage",
                    "cpswitch_c_active": False,
                }
                self.add_section(
                    name=f"XY coverage{title_suf}",
                    anchor=f"mosdepth-xy-coverage{id_suf}",
                    plot=bargraph.plot(xy_cov, xy_keys, pconfig),
                )

            if cumcov_dist_data:
                threshs, hidden_threshs = get_cov_thresholds()
                self.genstats_cov_thresholds(genstats, genstats_headers, cumcov_dist_data, threshs, hidden_threshs)
                self.genstats_mediancov(genstats, genstats_headers, cumcov_dist_data)

        # Add mean coverage to General Stats
        genstats_headers["mean_coverage"] = OrderedDict(
            {
                "title": "Mean Cov.",
                "description": "Mean coverage",
                "min": 0,
                "suffix": "X",
                "scale": "BuPu",
            }
        )
        self.general_stats_addcols(genstats, genstats_headers)

    def parse_cov_dist(self, scope):
        cumcov_dist_data = defaultdict(OrderedDict)  # cumulative distribution
        cov_dist_data = defaultdict(OrderedDict)  # absolute (non-cumulative) coverage
        perchrom_avg_data = defaultdict(OrderedDict)  # per chromosome average coverage

        # Parse coverage distributions
        for f in self.find_log_files(f"mosdepth/{scope}_dist"):
            s_name = self.clean_s_name(f["fn"], f)
            if s_name in cumcov_dist_data:  # both region and global might exist, prioritizing region
                continue

            for line in f["f"].split("\n"):
                if "\t" not in line:
                    continue
                contig, cutoff_reads, bases_fraction = line.split("\t")
                if float(bases_fraction) == 0:
                    continue

                # Parse cumulative coverage
                if contig == "total":
                    cumcov = 100.0 * float(bases_fraction)
                    x = int(cutoff_reads)
                    cumcov_dist_data[s_name][x] = cumcov

                # Calculate per-contig coverage
                else:
                    # filter out contigs based on exclusion patterns
                    if any(fnmatch.fnmatch(contig, str(pattern)) for pattern in self.cfg["exclude_contigs"]):
                        try:
                            if self.cfg.get("show_excluded_debug_logs") is True:
                                log.debug(f"Skipping excluded contig '{contig}'")
                        except (AttributeError, KeyError):
                            pass
                        continue

                    # filter out contigs based on inclusion patterns
                    if len(self.cfg["include_contigs"]) > 0 and not any(
                        fnmatch.fnmatch(contig, pattern) for pattern in self.cfg["include_contigs"]
                    ):
                        # Commented out since this could be many thousands of contigs!
                        # log.debug(f"Skipping not included contig '{contig}'")
                        continue

                    avg = perchrom_avg_data[s_name].get(contig, 0) + float(bases_fraction)
                    perchrom_avg_data[s_name][contig] = avg

            if s_name in cumcov_dist_data:
                self.add_data_source(f, s_name=s_name, section="genome_results")

        # Applying the contig coverage cutoff. First, count the total coverage for
        # every contig.
        total_cov_per_contig = defaultdict(lambda: 0)
        total_cov = 0
        for s_name, perchrom in perchrom_avg_data.items():
            for contig, cov in perchrom.items():
                total_cov_per_contig[contig] += cov
                total_cov += cov

        # Now, collecting the contigs that passed the cutoff.
        req_cov = float(total_cov) * self.cfg["perchrom_fraction_cutoff"]
        passing_contigs = set()
        for s_name, perchrom in perchrom_avg_data.items():
            for contig, cov in perchrom.items():
                if float(total_cov_per_contig[contig]) > req_cov:
                    if contig not in passing_contigs:
                        passing_contigs.add(contig)

        rejected_contigs = set()
        filtered_perchrom_avg_data = defaultdict(OrderedDict)
        for s_name, perchrom in perchrom_avg_data.items():
            for contig, cov in perchrom.items():
                if contig not in passing_contigs:
                    rejected_contigs.add(contig)
                else:
                    filtered_perchrom_avg_data[s_name][contig] = cov
        perchrom_avg_data = filtered_perchrom_avg_data

        if rejected_contigs:
            if self.cfg.get("show_excluded_debug_logs") is True:
                log.debug(
                    f"Skipping {len(rejected_contigs)} contigs not passing the "
                    f"cutoff of {self.cfg['perchrom_fraction_cutoff']}% of "
                    f"{total_cov:.2f}x total coverage, which is {req_cov:.2f}x. "
                    f"Skipping contigs: {''.join(rejected_contigs)}"
                )

        # Additionally, collect X and Y counts if we have them
        xy_cov = dict()
        for s_name, perchrom in perchrom_avg_data.items():
            x_cov = False
            y_cov = False
            for contig, cov in perchrom.items():
                if self.cfg.get("xchr"):
                    if str(self.cfg["xchr"]) == str(contig):
                        x_cov = cov
                else:
                    if contig.lower() == "x" or contig.lower() == "chrx":
                        x_cov = cov
                if self.cfg.get("ychr"):
                    if str(self.cfg["ychr"]) == str(contig):
                        y_cov = cov
                else:
                    if contig.lower() == "y" or contig.lower() == "chry":
                        y_cov = cov
            # Only save these counts if we have both x and y
            if x_cov and y_cov:
                xy_cov[s_name] = {"x": x_cov, "y": y_cov}

        # Correct per-contig average, since mosdepth reports cumulative coverage for at least
        # a certain value (see https://github.com/brentp/mosdepth#distribution-output).
        # For that reason, the 0 category (which is always 1) should not be included.
        for i in perchrom_avg_data:
            for j in perchrom_avg_data[i]:
                perchrom_avg_data[i][j] -= 1

        # Calculate absolute coverage distribution (global)
        for s_name, s_cumcov_dist in cumcov_dist_data.items():
            # Create sorted list of tuples (x, cumcov)
            cumcov_dist = sorted(s_cumcov_dist.items())

            # Calculate absolute coverage for the given x by taking the difference between
            # the current and previous cumulative coverage.
            #
            #   *example*              x:  cumcov:  abscov:
            #   3x                     3x  0      =               0
            #   2x     -               2x  0.10   = 0.10 - 0    = 0.10
            #   1x     --------        1x  0.80   = 0.80 - 0.10 = 0.70
            #   genome ..........      0x  1.00   = 1.00 - 0.80 = 0.20
            prev_x, prev_cumcov = cumcov_dist.pop()
            if not cumcov_dist:
                cov_dist_data[s_name][prev_x] = 1.0
            else:
                while cumcov_dist:
                    x, cumcov = cumcov_dist.pop()
                    cov_dist_data[s_name][x] = cumcov - prev_cumcov
                    prev_x, prev_cumcov = x, cumcov

        return cumcov_dist_data, cov_dist_data, perchrom_avg_data, xy_cov

    def genstats_cov_thresholds(self, genstats, genstats_headers, cumcov_dist_data, threshs, hidden_threshs):
        for s_name, d in cumcov_dist_data.items():
            dist_subset = {t: data for t, data in d.items() if t in threshs}
            for t in threshs:
                genstats[s_name][f"{t}_x_pc"] = dist_subset.get(t, 0)

        for t in threshs:
            genstats_headers[f"{t}_x_pc"] = {
                "title": f"&ge; {t}X",
                "description": f"Fraction of genome with at least {t}X coverage",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": t in hidden_threshs,
            }

    def genstats_mediancov(self, genstats, genstats_headers, cumcov_dist_data):
        for s_name, d in cumcov_dist_data.items():
            median_cov = None
            for this_cov, cum_pct in d.items():
                if cum_pct >= 50:
                    median_cov = this_cov
                    break
            genstats[s_name]["median_coverage"] = median_cov

        genstats_headers["median_coverage"] = {
            "title": "Median",
            "description": "Median coverage",
            "min": 0,
            "suffix": "X",
            "scale": "BuPu",
        }


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
