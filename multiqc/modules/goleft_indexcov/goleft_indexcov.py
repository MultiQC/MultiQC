"""MultiQC module to plot output from goleft indexcov

https://github.com/brentp/goleft/tree/master/indexcov
"""


import collections
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, scatter

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="goleft indexcov",
            anchor="goleft_indexcov",
            href="https://github.com/brentp/goleft/tree/master/indexcov",
            info="quickly estimates coverage from a whole-genome bam index.",
            doi="10.1093/gigascience/gix090",
        )

        # Parse ROC data
        self.roc_plot_data = collections.defaultdict(lambda: collections.defaultdict(dict))
        for f in self.find_log_files("goleft_indexcov/roc", filehandles=True):
            self.parse_roc_plot_data(f)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        num_roc_samples = 0
        for chrom in self.roc_plot_data:
            self.roc_plot_data[chrom] = self.ignore_samples(self.roc_plot_data[chrom])
            num_roc_samples = max(len(self.roc_plot_data[chrom]), num_roc_samples)

        # Write data to file
        self.write_data_file(self.roc_plot_data, "goleft_roc")

        # Parse BIN data
        self.bin_plot_data = {}
        for f in self.find_log_files("goleft_indexcov/ped", filehandles=True):
            self.parse_bin_plot_data(f)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.bin_plot_data = self.ignore_samples(self.bin_plot_data)

        # Write data to file
        self.write_data_file(self.roc_plot_data, "goleft_bin")

        # Stop execution if no samples
        num_samples = max(len(self.bin_plot_data), num_roc_samples)
        if num_samples == 0:
            raise UserWarning

        log.info("Found {} samples".format(num_samples))

        # Add sections to the report
        self.roc_plot()
        self.bin_plot()

    def _short_chrom(self, chrom):
        """Plot standard chromosomes + X, sorted numerically.

        Allows specification from a list of chromosomes via config
        for non-standard genomes.
        """
        default_allowed = set(["X"])
        allowed_chroms = set(getattr(config, "goleft_indexcov_config", {}).get("chromosomes", []))

        chrom_clean = chrom.replace("chr", "")
        try:
            chrom_clean = int(chrom_clean)
        except ValueError:
            if chrom_clean not in default_allowed and chrom_clean not in allowed_chroms:
                chrom_clean = None

        if allowed_chroms:
            if chrom in allowed_chroms or chrom_clean in allowed_chroms:
                return chrom_clean
        elif isinstance(chrom_clean, int) or chrom_clean in default_allowed:
            return chrom_clean

    def parse_roc_plot_data(self, f):
        header = f["f"].readline()
        sample_names = [self.clean_s_name(x, f) for x in header.strip().split()[2:]]
        for parts in (l.rstrip().split() for l in f["f"]):
            if len(parts) > 2:
                chrom, cov = parts[:2]
                sample_vals = parts[2:]
                if self._short_chrom(chrom) is not None:
                    for val, sample in zip(sample_vals, sample_names):
                        if not self.is_ignore_sample(sample):
                            self.roc_plot_data[chrom][sample][float(cov)] = float(val)

    def roc_plot(self):
        def to_padded_str(x):
            x = self._short_chrom(x)
            try:
                return "%06d" % x
            except TypeError:
                return x

        chroms = sorted(self.roc_plot_data.keys(), key=to_padded_str)

        max_chroms = getattr(config, "goleft_indexcov_config", {}).get("max_chroms", 50)
        if len(chroms) > max_chroms:
            log.warning("Too many chromosomes found: %s, limiting to %s" % (len(chroms), max_chroms))
            chroms = chroms[:max_chroms]

        pconfig = {
            "id": "goleft_indexcov-roc-plot",
            "title": "goleft indexcov: ROC - genome coverage per scaled depth by chromosome",
            "xlab": "Scaled coverage",
            "ylab": "Proportion of regions covered",
            "ymin": 0,
            "ymax": 1.0,
            "xmin": 0,
            "xmax": 1.5,
            "data_labels": [{"name": self._short_chrom(c)} for c in chroms],
        }
        self.add_section(
            name="Scaled coverage ROC plot",
            anchor="goleft_indexcov-roc",
            description="Coverage (ROC) plot that shows genome coverage at at given (scaled) depth.",
            helptext="""
                Lower coverage samples have shorter curves where the proportion of regions covered
                drops off more quickly. This indicates a higher fraction of low coverage regions.
            """,
            plot=linegraph.plot([self.roc_plot_data[c] for c in chroms], pconfig),
        )

    def parse_bin_plot_data(self, f):
        header = f["f"].readline()[1:].strip().split("\t")
        self.bin_plot_data_empty_samples = []
        for sample_parts in (l.split("\t") for l in f["f"]):
            cur = dict(zip(header, sample_parts))
            cur["sample_id"] = self.clean_s_name(cur["sample_id"], f)
            total = float(cur["bins.in"]) + float(cur["bins.out"])
            try:
                self.bin_plot_data[cur["sample_id"]] = {
                    "x": float(cur["bins.lo"]) / total,
                    "y": float(cur["bins.out"]) / total,
                }
            except ZeroDivisionError:
                self.bin_plot_data_empty_samples.append(cur["sample_id"])

    def bin_plot(self):
        pconfig = {
            "id": "goleft_indexcov-bin-plot",
            "title": "goleft indexcov: Problematic low and non-uniform coverage bins",
            "xlab": "Proportion of bins with depth < 0.15",
            "ylab": "Proportion of bins with depth outside of (0.85, 1.15)",
            "yCeiling": 1.0,
            "yFloor": 0.0,
            "xCeiling": 1.0,
            "xFloor": 0.0,
        }
        extra = ""
        if len(self.bin_plot_data_empty_samples) > 0:
            # Bootstrap alert about missing samples
            extra = f"""<div class="alert alert-warning" style="margin:2rem 0;">
                <strong>Warning:</strong>
                {len(self.bin_plot_data_empty_samples)} sample{'s' if len(self.bin_plot_data_empty_samples) > 1 else ''} had zero bins and could not be plotted.
                <a href="#goleft_empty_samples" onclick="$('#goleft_empty_samples').slideToggle();">Click to show missing sample names.</a>
                <div id="goleft_empty_samples" style="display:none;">
                    <ul><li><code>{"</code></li>, <li><code>".join(self.bin_plot_data_empty_samples)}</code></li></ul>
                </div>
            </div>"""
        self.add_section(
            name="Problem coverage bins",
            anchor="goleft_indexcov-bin",
            description="This plot identifies problematic samples using binned coverage distributions.",
            helptext="""
                We expect bins to be around 1, so deviations from this indicate problems.
                Low coverage bins (`< 0.15`) on the x-axis have regions with low or missing coverage.
                Higher values indicate truncated BAM files or missing data.
                Bins with skewed distributions (`<0.85` or `>1.15`) on the y-axis detect dosage bias.
                Large values on the y-axis are likely to impact CNV and structural variant calling.
                See the [goleft indexcov bin documentation](https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md)
                for more details.
            """,
            plot=scatter.plot(self.bin_plot_data, pconfig),
            content=extra,
        )
