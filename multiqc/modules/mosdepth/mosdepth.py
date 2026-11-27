import fnmatch
import gzip
import logging
import os
import zlib
from collections import defaultdict
from collections.abc import Iterable
from typing import cast

from multiqc import Plot, config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.qualimap.QM_BamQC import genome_fraction_helptext
from multiqc.plots import bargraph, linegraph
from multiqc.plots.linegraph import smooth_array
from multiqc.plots.table_object import ValueT
from multiqc.utils.util_functions import update_dict

log = logging.getLogger(__name__)


def read_config():
    cfg = getattr(config, "mosdepth_config", {})
    if not isinstance(cfg, dict):
        return {}

    cfg["include_contigs"] = cfg.get("include_contigs", [])
    if not isinstance(cfg["include_contigs"], list):
        cfg["include_contigs"] = []

    cfg["exclude_contigs"] = cfg.get("exclude_contigs", [])
    if not isinstance(cfg["exclude_contigs"], list):
        cfg["exclude_contigs"] = []

    xchr = cfg.get("xchr")
    if xchr and isinstance(cfg["xchr"], str):
        cfg["xchr"] = xchr

    ychr = cfg.get("ychr")
    if ychr and isinstance(cfg["ychr"], str):
        cfg["ychr"] = ychr

    if cfg["include_contigs"]:
        log.debug(f"Trying to include these contigs in mosdepth: {', '.join(cfg['include_contigs'])}")
    if cfg["exclude_contigs"]:
        log.debug(f"Excluding these contigs from mosdepth: {', '.join(cfg['exclude_contigs'])}")
    if cfg.get("xchr"):
        log.debug(f'Using "{cfg["xchr"]}" as X chromosome name')
    if cfg.get("ychr"):
        log.debug(f'Using "{cfg["ychr"]}" as Y chromosome name')

    cutoff = cfg.get("perchrom_fraction_cutoff", 0.0)
    try:
        cutoff = float(cutoff)
    except ValueError:
        cutoff = 0.0
    if cutoff != 0.0:
        log.debug(f"Setting mosdepth coverage cutoff to display the contigs to {cutoff * 100.0}%")
    cfg["perchrom_fraction_cutoff"] = cutoff

    pass_pct = cfg.get("region_pass_pct", DEFAULT_REGION_PASS_PCT)
    try:
        pass_pct = float(pass_pct)
    except (TypeError, ValueError):
        pass_pct = -1.0
    if not 0 < pass_pct <= 100:
        log.warning(
            f"mosdepth_config.region_pass_pct must be a percentage above 0 and at most 100, "
            f"got {cfg['region_pass_pct']!r}. Using {DEFAULT_REGION_PASS_PCT}."
        )
        pass_pct = DEFAULT_REGION_PASS_PCT
    cfg["region_pass_pct"] = pass_pct

    return cfg


DEFAULT_REGION_PASS_PCT = 99.0

# Errors from reading a file that matches *.thresholds.bed.gz but isn't usable mosdepth output:
# bad content, not gzip, or truncated (UnicodeDecodeError is a ValueError, gzip.BadGzipFile an OSError)
THRESHOLDS_PARSE_ERRORS = (ValueError, OSError, EOFError, zlib.error)


def count_passing_regions(lines: Iterable[str], pass_pct: float) -> tuple[int, dict[int, int]]:
    """
    Count, for each threshold, the regions with at least `pass_pct` % of bases covered at that
    threshold, from the lines of a decompressed {prefix}.thresholds.bed.gz:

    #chrom  start  end  region  1X  10X  20X  30X  (threshold columns follow --thresholds)

    Streams over the lines and keeps only counters, so memory does not grow with the number of
    regions. Zero-length regions are left out. Returns the number of regions counted and the
    number passing per threshold. Raises ValueError if the content is not mosdepth thresholds output.
    """
    lines_iter = iter(lines)
    header_line = next(lines_iter, None)
    if header_line is None:
        raise ValueError("empty file")
    header = header_line.rstrip("\r\n").lstrip("#").split("\t")
    if header[:4] != ["chrom", "start", "end", "region"] or not all(c.endswith("X") for c in header[4:]):
        raise ValueError(f"unexpected header: {header_line.strip()!r}")
    thresholds = [int(c[:-1]) for c in header[4:]]

    n_regions = 0
    passing = dict.fromkeys(thresholds, 0)
    for line_num, line in enumerate(lines_iter, start=2):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != len(header):
            raise ValueError(f"line {line_num} has {len(fields)} columns, expected {len(header)} from the header")
        length = int(fields[2]) - int(fields[1])
        if length <= 0:
            continue
        n_regions += 1
        for t, count in zip(thresholds, fields[4:]):
            # Same as count / length * 100 >= pass_pct, without dividing
            if int(count) * 100 >= pass_pct * length:
                passing[t] += 1
    return n_regions, passing


def genstats_cov_thresholds(cum_fraction_by_cov: dict[int, float], threshs: list[int]) -> dict[str, float]:
    genstats: dict[str, float] = {}
    sorted_cum_fraction_by_cov = sorted(cum_fraction_by_cov.items())
    for t in threshs:
        # take next known value
        # e.g. if only 50x is known but the threshold is 40x, take the 50x value
        # if there is no next threshold, just take 0.0
        cov_val = next((proportion for cov, proportion in sorted_cum_fraction_by_cov if cov >= t), 0.0)
        genstats[f"{t}_x_pc"] = cov_val * 100.0
    return genstats


def cov_at_cum_fraction(cum_fraction_by_cov: dict[int, float], fraction: float) -> int | None:
    """Highest coverage at which at least `fraction` of bases are still covered."""
    for cov, cum_fraction in sorted(cum_fraction_by_cov.items(), reverse=True):
        if cum_fraction >= fraction:
            return cov
    return None


def calc_iqr_coverage(cum_fraction_by_cov: dict[int, float]) -> int | None:
    """Interquartile range of coverage.

    The fractions are cumulative ("at least X coverage"), so Q3 is where 25% of
    bases remain covered and Q1 is where 75% remain.
    """
    q3 = cov_at_cum_fraction(cum_fraction_by_cov, 0.25)
    q1 = cov_at_cum_fraction(cum_fraction_by_cov, 0.75)
    return None if q1 is None or q3 is None else q3 - q1


class MultiqcModule(BaseMultiqcModule):
    """
    Mosdepth can generate several output files all with a common prefix and different endings:

    - per-base depth (`{prefix}.per-base.bed.gz`),
    - mean per-window depth given a window size (`{prefix}.regions.bed.gz`, if a BED file provided with `--by`),
    - mean per-region given a BED file of regions (`{prefix}.regions.bed.gz`, if a window size provided with `--by`),
    - a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide (`{prefix}.mosdepth.global.dist.txt` and `{prefix}.mosdepth.region.dist.txt`),
    - quantized output that merges adjacent bases as long as they fall in the same coverage bins (`{prefix}.quantized.bed.gz`),
    - threshold output to indicate how many bases in each region are covered at the given thresholds (`{prefix}.thresholds.bed.gz`)
    - summary output providing region length, coverage mean, min, and max for each region. (`{prefix}.mosdepth.summary.txt`)

    The MultiQC module plots coverage distributions from 2 kinds of outputs:

    - `{prefix}.mosdepth.region.dist.txt`
    - `{prefix}.mosdepth.global.dist.txt`

    Using "region" if exists, otherwise "global". Plotting 3 figures:

    - Proportion of bases in the reference genome with, at least, a given depth of coverage (cumulative coverage distribution).
    - Proportion of bases in the reference genome with a given depth of coverage (absolute coverage distribution).
    - Average coverage per contig/chromosome.

    Also plotting the percentage of the genome covered at a threshold in the General Stats section.
    The default thresholds are 1, 5, 10, 30, 50, which can be customised in the config as follows:

    ```yaml
    mosdepth_config:
      general_stats_coverage:
        - 10
        - 20
        - 40
        - 200
        - 30000
    ```

    You can also specify which columns would be hidden when the report loads (by default, all values are hidden except 30X):

    ```yaml
    general_stats_coverage_hidden:
      - 10
      - 20
      - 200
    ```

    A hidden `IQR/Median` General Stats column measures coverage uniformity: the interquartile
    range of the coverage distribution divided by the median coverage. Lower values mean more
    uniform coverage.

    If mosdepth was run with `--by` and `--thresholds`, the `{prefix}.thresholds.bed.gz` file is
    used to add a pair of General Stats columns per threshold X: the % of regions with at least
    99% of their bases covered at X or more, and the number of regions that miss this. The
    `% bases ≥ X` columns are pooled over all target bases, so they look the same whether the
    poorly covered bases are spread thinly over many regions or concentrated in one region with
    no coverage at all; the failing-region count tells those two cases apart. Columns follow the
    same visibility as the `% bases ≥ X` columns: only thresholds that are in `general_stats_coverage`
    and not in `general_stats_coverage_hidden` are shown by default. The required % of bases can
    be changed in the config:

    ```yaml
    mosdepth_config:
      # Percentage of a region's bases that must reach the threshold for the region to pass
      region_pass_pct: 95
    ```

    Note that `*.thresholds.bed.gz` is matched on the filename only, since MultiQC can't read
    the contents of gzipped files during the file search. A file with this name that isn't
    mosdepth output is skipped with a warning.

    For the per-contig coverage plot, you can include and exclude contigs based on name or pattern.

    For example, you could add the following to your MultiQC config file:

    ```yaml
    mosdepth_config:
      include_contigs:
        - "chr*"
      exclude_contigs:
        - "*_alt"
        - "*_decoy"
        - "*_random"
        - "chrUn*"
        - "HLA*"
        - "chrM"
        - "chrEBV"
    ```

    Note that exclusion superseeds inclusion for the contig filters.

    To additionally avoid cluttering the plot, mosdepth can exclude contigs with a low relative coverage.

    ```yaml
    mosdepth_config:
      # Should be a fraction, e.g. 0.001 (exclude contigs with 0.1% coverage of sum of
      # coverages across all contigs)
      perchrom_fraction_cutoff: 0.001
    ```

    If you want to see what is being excluded, you can set `show_excluded_debug_logs` to `True`:

    ```yaml
    mosdepth_config:
      show_excluded_debug_logs: True
    ```

    This will then print a debug log message (use `multiqc -v`) for each excluded contig.
    This is disabled by default as there can be very many in some cases.

    Besides the `{prefix}.mosdepth.global.dist.txt` and `{prefix}.mosdepth.region.dist.txt`
    files, the `{prefix}.mosdepth.summary.txt` file is used for the General Stats table.

    The module also plots an X/Y relative chromosome coverage per sample. By default, it finds chromosome named X/Y or chrX/chrY, but that can be customised:

    ```yaml
    mosdepth_config:
      # Name of the X and Y chromosomes. If not specified, MultiQC will search for
      # any chromosome names that look like x, y, chrx or chry (case-insensitive)
      xchr: myXchr
      ychr: myYchr
    ```
    """

    def __init__(self):
        super().__init__(
            name="Mosdepth",
            anchor="mosdepth",
            href="https://github.com/brentp/mosdepth",
            info="Fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing",
            doi="10.1093/bioinformatics/btx699",
            license="MIT License",
            license_url="https://github.com/brentp/mosdepth/blob/master/LICENSE",
        )

        self.cfg = read_config()
        genstats_by_sample: dict[str, dict[str, int | float]] = defaultdict(dict)  # mean coverage

        # Parse mean coverage
        for f in self.find_log_files("mosdepth/summary"):
            s_name = self.clean_s_name(f["fn"], f)
            sample_stats: dict[str, int | float] = {}
            try:
                for line in f["f"].splitlines():
                    # The first column can be a contig name, "total", "total_region".
                    # We want to use "total_region" if available. It is available when
                    # --by is specified. It always goes after "total", so we can just
                    # assume it will override the information collected for "total":
                    if line.startswith(("total\t", "total_region\t")):
                        _contig, length, bases, mean, min_cov, max_cov = line.split("\t")
                        sample_stats["mean_coverage"] = float(mean)
                        sample_stats["min_coverage"] = float(min_cov)
                        sample_stats["max_coverage"] = float(max_cov)
                        sample_stats["coverage_bases"] = int(bases)
                        sample_stats["length"] = int(length)
            except ValueError:
                log.warning(
                    f"Skipping {f['root']}/{f['fn']}: could not parse line {line!r}. "
                    f"The file may be truncated or corrupted."
                )
                continue
            if sample_stats:
                genstats_by_sample[s_name].update(sample_stats)
                self.add_data_source(f, s_name=s_name, section="summary")

        # Filter out any samples from --ignore-samples
        genstats_by_sample = defaultdict(dict, self.ignore_samples(genstats_by_sample))
        samples_in_summary = set(genstats_by_sample.keys())

        raw_cov_dist_global = self.parse_cov_dist("global")
        raw_cov_dist_region = self.parse_cov_dist("region")
        data_dicts_global = [self.ignore_samples(d) for d in raw_cov_dist_global]
        data_dicts_region = [self.ignore_samples(d) for d in raw_cov_dist_region]

        samples_global = set.union(*(set(d.keys()) for d in data_dicts_global))
        samples_region = set.union(*(set(d.keys()) for d in data_dicts_region))
        regions_by_sample, region_threshs = self.parse_thresholds_bed()
        regions_by_sample = self.ignore_samples(regions_by_sample)
        update_dict(genstats_by_sample, regions_by_sample)

        samples_found = samples_in_summary | samples_global | samples_region | set(regions_by_sample)
        if not samples_found:
            raise ModuleNoSamplesFound
        log.info(f"Found reports for {len(samples_found)} samples")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        descr_suf = ""
        if any(data_dicts_global) and any(data_dicts_region):
            descr_suf = (
                f". Note that for {len(samples_global)} samples, a BED file was "
                f"provided, so the data was calculated across those regions. For "
                f"{len(samples_region)} samples, it's calculated across "
                f"the entire genome length"
            )
            samples_both = samples_global & samples_region
            if samples_both:
                descr_suf += (
                    f". {len(samples_both)} samples have both global and region "
                    f"reports, and we are showing the data for regions"
                )
        elif samples_region:
            descr_suf = ". Calculated across the target regions"
        elif samples_global:
            descr_suf = ". Calculated across the entire genome length"

        if samples_region or samples_global:
            # Prioritizing region reports if found
            data_dicts = data_dicts_global
            for d, d_region in zip(data_dicts, data_dicts_region):
                d.update(d_region)
            (
                cum_cov_dist_by_sample,
                perchrom_avg_by_sample,
                xy_cov_by_sample,
                extra_genstats_by_sample,
            ) = data_dicts

            if cum_cov_dist_by_sample:
                xmax = 0
                for cum_cov_by_x in cum_cov_dist_by_sample.values():
                    for x, cumcov in cum_cov_by_x.items():
                        if cumcov is not None and cumcov > 1:  # require >1% to prevent long flat tail
                            xmax = max(xmax, x)

                # Write data to file, sort columns numerically and convert to strings
                cumcov_dist_data_writeable = {
                    sample: {str(k): v for k, v in sorted(cum_cov_by_x.items())}
                    for sample, cum_cov_by_x in cum_cov_dist_by_sample.items()
                }
                self.write_data_file(cumcov_dist_data_writeable, "mosdepth_cumcov_dist")

                self.add_section(
                    name="Cumulative coverage distribution",
                    anchor="mosdepth-cumcoverage-dist",
                    description=(
                        f"Proportion of bases in the reference genome with, "
                        f"at least, a given depth of coverage{descr_suf}"
                    ),
                    helptext=genome_fraction_helptext,
                    plot=linegraph.plot(
                        cum_cov_dist_by_sample,
                        {
                            "id": "mosdepth-cumcoverage-dist-id",
                            "title": "Mosdepth: Cumulative coverage distribution",
                            "xlab": "Cumulative Coverage (X)",
                            "ylab": "% bases in genome/regions covered by at least X reads",
                            "ymin": 0,
                            "ymax": 100,
                            "xmin": 0,
                            "xmax": xmax,
                            "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
                            "smooth_points": 500,
                        },
                    ),
                )

                # Write data to file, sort columns numerically and convert to strings
                cov_dist_data_writeable = {
                    sample: {str(k): v for k, v in sorted(cum_cov_by_x.items())}
                    for sample, cum_cov_by_x in cum_cov_dist_by_sample.items()
                }
                self.write_data_file(cov_dist_data_writeable, "mosdepth_cov_dist")

            if perchrom_avg_by_sample:
                # Write data to file
                self.write_data_file(perchrom_avg_by_sample, "mosdepth_perchrom")

                num_contigs = max([len(x.keys()) for x in perchrom_avg_by_sample.values()])
                perchrom_plot: Plot | str | None
                if num_contigs > 1:
                    perchrom_plot = linegraph.plot(
                        perchrom_avg_by_sample,
                        {
                            "id": "mosdepth-coverage-per-contig-multi",
                            "title": "Mosdepth: Coverage per contig",
                            "xlab": "Region",
                            "ylab": "Average Coverage",
                            "tt_decimals": 1,
                            "tt_suffix": "x",
                            "smooth_points": 500,
                            "logswitch": True,
                            "categories": True,
                        },
                    )
                else:
                    perchrom_plot = bargraph.plot(
                        perchrom_avg_by_sample,
                        pconfig={
                            "id": "mosdepth-coverage-per-contig-single",
                            "title": "Mosdepth: Coverage per contig",
                            "xlab": "Sample",
                            "ylab": "Average Coverage",
                            "tt_suffix": "x",
                            "hide_zero_cats": False,
                            "categories": True,
                        },
                    )

                self.add_section(
                    name="Average coverage per contig",
                    anchor="mosdepth-coverage-per-contig-section",
                    description="Average coverage per contig or chromosome",
                    plot=perchrom_plot,
                )

            if xy_cov_by_sample:
                xy_keys = {
                    "x": {"name": self.cfg.get("xchr") or "Chromosome X"},
                    "y": {"name": self.cfg.get("ychr") or "Chromosome Y"},
                }
                pconfig = {
                    "id": "mosdepth-xy-coverage-plot",
                    "title": "Mosdepth: chrXY coverage",
                    "ylab": "Coverage",
                    "ysuffix": "X",
                    "cpswitch_c_active": False,
                }
                self.add_section(
                    name="XY coverage",
                    anchor="mosdepth-xy-coverage",
                    plot=bargraph.plot(
                        {sname: {"x": x_cov, "y": y_cov} for sname, (x_cov, y_cov) in xy_cov_by_sample.items()},
                        xy_keys,
                        pconfig,
                    ),
                )

            if extra_genstats_by_sample:
                update_dict(genstats_by_sample, extra_genstats_by_sample)

        genstats_headers = {}
        threshs, hidden_threshs = config.get_cov_thresholds("mosdepth_config")
        for t in threshs:
            genstats_headers[f"{t}_x_pc"] = {
                "title": f"≥ {t}X",
                "description": f"Fraction of genome with at least {t}X coverage",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": t in hidden_threshs,
            }
        pass_pct = self.cfg["region_pass_pct"]
        for t in sorted(region_threshs):
            hidden = t not in threshs or t in hidden_threshs
            genstats_headers[f"{t}_x_regions_pc"] = {
                "title": f"Regions ≥ {t}X",
                "description": f"% of regions with ≥{pass_pct:g}% of bases covered at ≥{t}X",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": hidden,
            }
            genstats_headers[f"{t}_x_regions_failed"] = {
                "title": f"Regions < {t}X",
                "description": f"Number of regions with less than {pass_pct:g}% of bases covered at ≥{t}X",
                "min": 0,
                "scale": "OrRd",
                "format": "{:,d}",
                "hidden": hidden,
            }
        # Add mosdepth summary to General Stats
        genstats_headers.update(
            {
                "median_coverage": {
                    "title": "Median",
                    "description": "Median coverage",
                    "min": 0,
                    "suffix": "X",
                    "scale": "BuPu",
                },
                "mean_coverage": {
                    "title": "Mean Cov.",
                    "description": "Mean coverage",
                    "min": 0,
                    "suffix": "X",
                    "scale": "BuPu",
                },
                "min_coverage": {
                    "title": "Min Cov.",
                    "description": "Minimum coverage",
                    "min": 0,
                    "suffix": "X",
                    "scale": "BuPu",
                    "hidden": True,
                },
                "max_coverage": {
                    "title": "Max Cov.",
                    "description": "Maximum coverage",
                    "min": 0,
                    "suffix": "X",
                    "scale": "BuPu",
                    "hidden": True,
                },
                "coverage_bases": {
                    "title": f"{config.base_count_prefix} Total Coverage Bases",
                    "description": f"Total coverage of bases ({config.base_count_desc})",
                    "min": 0,
                    "shared_key": "base_count",
                    "scale": "Greens",
                    "hidden": True,
                },
                "length": {
                    "title": "Genome length",
                    "description": "Total length of the genome",
                    "min": 0,
                    "scale": "Greys",
                    "format": "{:,d}",
                    "hidden": True,
                },
                "regions_total": {
                    "title": "Regions",
                    "description": "Number of regions in the thresholds.bed.gz file (zero-length regions excluded)",
                    "min": 0,
                    "scale": "Greys",
                    "format": "{:,d}",
                    "hidden": True,
                },
                "coverage_iqr_cv": {
                    "title": "IQR/Median",
                    "description": "Coverage IQR / median coverage",
                    "min": 0,
                    "scale": "BuPu",
                    "format": "{:,.2f}",
                    "hidden": True,
                },
            },
        )
        self.general_stats_addcols(cast(dict[str, dict[str, ValueT]], genstats_by_sample), genstats_headers)

    def parse_thresholds_bed(self) -> tuple[dict[str, dict[str, int | float]], set[int]]:
        """
        Summarise each {prefix}.thresholds.bed.gz (mosdepth `--by ... --thresholds ...`) into
        per-sample General Stats values: per threshold, the % of regions passing and the number
        of failing regions. Also returns the union of thresholds seen across samples.
        """
        pass_pct = self.cfg["region_pass_pct"]
        stats_by_sample: dict[str, dict[str, int | float]] = {}
        all_threshs: set[int] = set()
        for f in self.find_log_files("mosdepth/thresholds_bed", filecontents=False, filehandles=False):
            s_name = self.clean_s_name(f["fn"], f)
            path = os.path.join(f["root"], f["fn"])
            # The pattern can only match on the filename (file search can't look inside gzip),
            # so a same-named file from another tool is skipped rather than crashing the module.
            try:
                with gzip.open(path, "rt") as fh:
                    n_regions, passing = count_passing_regions(fh, pass_pct)
            except THRESHOLDS_PARSE_ERRORS as e:
                log.warning(f"Skipping {path}: not usable mosdepth thresholds.bed.gz output ({e})")
                continue
            if n_regions == 0:
                log.warning(f"Skipping {path}: no regions with non-zero length")
                continue

            stats: dict[str, int | float] = {"regions_total": n_regions}
            for t, n_pass in passing.items():
                # >=0X always passes, so the column carries no information
                if t == 0:
                    continue
                stats[f"{t}_x_regions_pc"] = 100.0 * n_pass / n_regions
                stats[f"{t}_x_regions_failed"] = n_regions - n_pass
                all_threshs.add(t)
            stats_by_sample[s_name] = stats
            self.add_data_source(f, s_name=s_name, section="thresholds")
        return stats_by_sample, all_threshs

    def parse_cov_dist(
        self, scope: str
    ) -> tuple[
        dict[str, dict],
        dict[str, dict],
        dict[str, tuple[float, float]],
        dict[str, dict[str, float | int | None]],
    ]:
        """
        Two types of coverage distributions are parsed: global and region.

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
        """

        cumulative_pct_by_cov_by_sample: dict[str, dict] = defaultdict(dict)  # cumulative distribution
        bases_fraction_sum_per_contig_per_sample: dict[str, dict[str, float]] = defaultdict(
            dict
        )  # per chromosome average coverage
        xy_cov_by_sample: dict[str, tuple[float, float]] = {}
        genstats_by_sample: dict[str, dict[str, float | int | None]] = {}

        threshs, _hidden_threshs = config.get_cov_thresholds("mosdepth_config")

        excluded_contigs = set()
        included_contigs = set()
        show_excluded_debug_logs = self.cfg.get("show_excluded_debug_logs") is True

        # Parse coverage distributions
        for f in self.find_log_files(f"mosdepth/{scope}_dist", filecontents=False, filehandles=True):
            s_name = self.clean_s_name(f["fn"], f)
            if s_name in cumulative_pct_by_cov_by_sample:  # both region and global might exist, prioritizing region
                continue

            bases_fraction_sum_per_contig: dict[str, float] = defaultdict(float)
            cum_fraction_by_cov: dict[int, float] = {}

            line_num = 0
            try:
                for line_num, line in enumerate(f["f"], start=1):
                    contig, cutoff_reads, bases_fraction = str(line).rstrip("\r\n").split("\t")
                    if bases_fraction == "0.00":
                        continue

                    # Parse cumulative coverage
                    if contig == "total":
                        cum_fraction_by_cov[int(cutoff_reads)] = float(bases_fraction)

                    # Calculate per-contig coverage
                    else:
                        # filter out contigs based on exclusion patterns
                        if contig in excluded_contigs:
                            continue

                        if contig not in included_contigs:
                            if any(fnmatch.fnmatch(contig, str(pattern)) for pattern in self.cfg["exclude_contigs"]):
                                excluded_contigs.add(contig)
                                if show_excluded_debug_logs:
                                    log.debug(f"Skipping excluded contig '{contig}'")
                                continue

                            # filter out contigs based on inclusion patterns
                            if len(self.cfg["include_contigs"]) > 0 and not any(
                                fnmatch.fnmatch(contig, pattern) for pattern in self.cfg["include_contigs"]
                            ):
                                # Commented out since this could be many thousands of contigs!
                                # log.debug(f"Skipping not included contig '{contig}'")
                                continue

                            included_contigs.add(contig)

                        bases_fraction_sum_per_contig[contig] += float(bases_fraction)
            except ValueError:
                log.warning(
                    f"Skipping {f['root']}/{f['fn']}: could not parse line {line_num}. The file may be truncated or corrupted."
                )
                continue
            if not cum_fraction_by_cov:
                log.warning(f"Skipping {f['root']}/{f['fn']}: no 'total' rows found. The file may be truncated.")
                continue

            self.add_data_source(f, s_name=s_name, section="genome_results")

            genstats_by_sample[s_name] = {}
            for k, v in genstats_cov_thresholds(cum_fraction_by_cov, threshs).items():
                genstats_by_sample[s_name][k] = v
            median_coverage = cov_at_cum_fraction(cum_fraction_by_cov, 0.5)
            genstats_by_sample[s_name]["median_coverage"] = median_coverage
            iqr_coverage = calc_iqr_coverage(cum_fraction_by_cov)
            if iqr_coverage is not None and median_coverage not in (None, 0):
                genstats_by_sample[s_name]["coverage_iqr_cv"] = iqr_coverage / median_coverage

            # Downsampling the data to avoid carrying a lot for the line plot that would downsample anyway
            cum_fraction_by_cov = dict(smooth_array(list(cum_fraction_by_cov.items()), 500))
            cumulative_pct_by_cov_by_sample[s_name] = {
                cutoff_reads: 100.0 * bases_fraction for cutoff_reads, bases_fraction in cum_fraction_by_cov.items()
            }
            bases_fraction_sum_per_contig_per_sample[s_name] = bases_fraction_sum_per_contig

        # Applying the contig coverage cutoff. First, count the total coverage for
        # every contig.
        total_cov_per_contig: dict[str, float] = defaultdict(lambda: 0)
        total_cov = 0.0
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                total_cov_per_contig[contig] += bases_fraction_sum
                total_cov += bases_fraction_sum

        # Now, collecting the contigs that passed the cutoff.
        req_cov = float(total_cov) * self.cfg["perchrom_fraction_cutoff"]
        passing_contigs = set()
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                if float(total_cov_per_contig[contig]) > req_cov and contig not in passing_contigs:
                    passing_contigs.add(contig)

        rejected_contigs = set()
        filtered_perchrom_avg_data: dict[str, dict[str, float]] = defaultdict(dict)
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                if contig not in passing_contigs:
                    rejected_contigs.add(contig)
                else:
                    filtered_perchrom_avg_data[s_name][contig] = bases_fraction_sum
        bases_fraction_sum_per_contig_per_sample = filtered_perchrom_avg_data

        if rejected_contigs and self.cfg.get("show_excluded_debug_logs") is True:
            log.debug(
                f"Skipping {len(rejected_contigs)} contigs not passing the "
                f"cutoff of {self.cfg['perchrom_fraction_cutoff']}% of "
                f"{total_cov:.2f}x total coverage, which is {req_cov:.2f}x. "
                f"Skipping contigs: {''.join(rejected_contigs)}"
            )

        # Additionally, collect X and Y counts if we have them
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            x_cov: float | None = None
            y_cov: float | None = None
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                if self.cfg.get("xchr"):
                    if str(self.cfg["xchr"]) == str(contig):
                        x_cov = bases_fraction_sum
                else:
                    if contig.lower() == "x" or contig.lower() == "chrx":
                        x_cov = bases_fraction_sum
                if self.cfg.get("ychr"):
                    if str(self.cfg["ychr"]) == str(contig):
                        y_cov = bases_fraction_sum
                else:
                    if contig.lower() == "y" or contig.lower() == "chry":
                        y_cov = bases_fraction_sum
            # Only save these counts if we have both x and y
            if x_cov and y_cov:
                xy_cov_by_sample[s_name] = x_cov, y_cov

        # Correct per-contig average, since mosdepth reports cumulative coverage for at least
        # a certain value (see https://github.com/brentp/mosdepth#distribution-output).
        # For that reason, the 0 category (which is always 1) should not be included.
        for s_name in bases_fraction_sum_per_contig_per_sample:
            for contig in bases_fraction_sum_per_contig_per_sample[s_name]:
                bases_fraction_sum_per_contig_per_sample[s_name][contig] -= 1

        return (
            cumulative_pct_by_cov_by_sample,
            bases_fraction_sum_per_contig_per_sample,
            xy_cov_by_sample,
            genstats_by_sample,
        )
