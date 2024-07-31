from collections import defaultdict
import logging
from typing import Optional, Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, scatter

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module uses the PED and ROC data files to create diagnostic plots of coverage per
    sample, helping to identify sample gender and coverage issues.

    By default, we attempt to only plot chromosomes using standard human-like naming
    (chr1, chr2... chrX or 1, 2 ... X) but you can specify chromosomes for detailed
    ROC plots for alternative naming schemes in your configuration with:

    ```yaml
    goleft_indexcov_config:
      chromosomes:
        - I
        - II
        - III
    ```

    The number of plotted chromosomes is limited to 50 by default, you can customise this with the following:

    ```yaml
    goleft_indexcov_config:
      max_chroms: 80
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="goleft indexcov",
            anchor="goleft_indexcov",
            href="https://github.com/brentp/goleft/tree/master/indexcov",
            info="Quickly estimate coverage from a whole-genome bam index, providing 16KB resolution",
            extra="This is useful as a quick QC to get coverage values across the genome.",
            doi="10.1093/gigascience/gix090",
        )

        # Parse ROC data
        roc_plot_data: Dict[str, Dict[str, Dict[float, float]]] = defaultdict(lambda: defaultdict(dict))
        for f in self.find_log_files("goleft_indexcov/roc", filehandles=True):
            header = f["f"].readline()
            sample_names = [self.clean_s_name(x, f) for x in header.strip().split()[2:]]
            one_file_data = defaultdict(lambda: defaultdict(dict))
            lines = f["f"].readlines()
            for line in lines:
                parts = line.rstrip().split()
                if len(parts) > 2:
                    chrom, cov = parts[:2]
                    sample_vals = parts[2:]
                    if self._short_chrom(chrom) is not None:
                        for val, sample in zip(sample_vals, sample_names):
                            if not self.is_ignore_sample(sample):
                                one_file_data[chrom][sample][float(cov)] = float(val)
            if not one_file_data:  # No chromosomes passed filter, attempt to parse all
                for line in lines:
                    parts = line.rstrip().split()
                    if len(parts) > 2:
                        chrom, cov = parts[:2]
                        sample_vals = parts[2:]
                        for val, sample in zip(sample_vals, sample_names):
                            if not self.is_ignore_sample(sample):
                                one_file_data[chrom][sample][cov] = val
            if one_file_data:
                self.add_data_source(f)
                for chrom, chrom_data in one_file_data.items():
                    for sn, sample_data in chrom_data.items():
                        for k, v in sample_data.items():
                            roc_plot_data[chrom][sn][k] = v

        # Filter to strip out ignored sample names
        num_roc_samples = 0
        for chrom in roc_plot_data:
            roc_plot_data[chrom] = self.ignore_samples(roc_plot_data[chrom])
            num_roc_samples = max(len(roc_plot_data[chrom]), num_roc_samples)

        # Write data to file
        self.write_data_file(roc_plot_data, "goleft_roc")

        # Parse BIN data
        bin_plot_data = {}
        bin_plot_data_empty_samples = []
        for f in self.find_log_files("goleft_indexcov/ped", filehandles=True):
            self.parse_bin_plot_data(f, bin_plot_data, bin_plot_data_empty_samples)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        bin_plot_data = self.ignore_samples(bin_plot_data)

        # Write data to file
        self.write_data_file(roc_plot_data, "goleft_bin")

        # Stop execution if no samples
        num_samples = max(len(bin_plot_data), num_roc_samples)
        if num_samples == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {num_samples} samples")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add sections to the report
        if roc_plot_data:
            self.roc_plot(roc_plot_data)
        if bin_plot_data:
            self.bin_plot(bin_plot_data, bin_plot_data_empty_samples)

    @staticmethod
    def _short_chrom(chrom: str) -> Optional[str]:
        """Plot standard chromosomes + X, sorted numerically.

        Allows specification from a list of chromosomes via config
        for non-standard genomes.
        """
        default_allowed = {"X"}
        allowed_chroms = [str(c) for c in set(getattr(config, "goleft_indexcov_config", {}).get("chromosomes", []))]

        chrom_clean = chrom.replace("chr", "")
        try:
            chrom_clean = int(chrom_clean)
        except ValueError:
            if chrom_clean not in default_allowed and chrom_clean not in allowed_chroms:
                chrom_clean = None
        else:
            chrom_clean = str(chrom_clean)

        if (
            allowed_chroms
            and chrom in allowed_chroms
            or chrom_clean in allowed_chroms
            or chrom_clean in default_allowed
        ):
            return chrom_clean
        return None

    def roc_plot(self, roc_plot_data):
        def to_padded_str(x):
            short_chrom = self._short_chrom(x)
            try:
                return "%06d" % short_chrom
            except TypeError:
                return x

        chroms = sorted(roc_plot_data.keys(), key=to_padded_str)

        max_chroms = getattr(config, "goleft_indexcov_config", {}).get("max_chroms", 50)
        if len(chroms) > max_chroms:
            log.warning(f"Too many chromosomes found: {len(chroms)}, limiting to {max_chroms}")
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
            "data_labels": [{"name": self._short_chrom(c) if self._short_chrom(c) else c} for c in chroms],
        }
        self.add_section(
            name="Scaled coverage ROC plot",
            anchor="goleft_indexcov-roc",
            description="Coverage (ROC) plot that shows genome coverage at given (scaled) depth.",
            helptext="""
                Lower coverage samples have shorter curves where the proportion of regions covered
                drops off more quickly. This indicates a higher fraction of low coverage regions.
            """,
            plot=linegraph.plot([roc_plot_data[c] for c in chroms], pconfig),
        )

    def parse_bin_plot_data(self, f, bin_plot_data, bin_plot_data_empty_samples):
        header = f["f"].readline()[1:].strip().split("\t")
        for sample_parts in (line.split("\t") for line in f["f"]):
            cur = dict(zip(header, sample_parts))
            cur["sample_id"] = self.clean_s_name(cur["sample_id"], f)
            total = float(cur["bins.in"]) + float(cur["bins.out"])
            try:
                bin_plot_data[cur["sample_id"]] = {
                    "x": float(cur["bins.lo"]) / total,
                    "y": float(cur["bins.out"]) / total,
                }
            except ZeroDivisionError:
                bin_plot_data_empty_samples.append(cur["sample_id"])

    def bin_plot(self, bin_plot_data, bin_plot_data_empty_samples):
        pconfig = {
            "id": "goleft_indexcov-bin-plot",
            "title": "goleft indexcov: Problematic low and non-uniform coverage bins",
            "xlab": "Proportion of bins with depth < 0.15",
            "ylab": "Proportion of bins with depth outside of (0.85, 1.15)",
            "ymax": 1.0,
            "ymin": 0.0,
            "x_clipmax": 1.0,
            "x_clipmin": 0.0,
        }
        extra = ""
        if len(bin_plot_data_empty_samples) > 0:
            # Bootstrap alert about missing samples
            extra = f"""<div class="alert alert-warning" style="margin:2rem 0;">
                <strong>Warning:</strong>
                {len(bin_plot_data_empty_samples)} sample{'s' if len(bin_plot_data_empty_samples) > 1 else ''} had zero bins and could not be plotted.
                <a href="#goleft_empty_samples" onclick="$('#goleft_empty_samples').slideToggle();">Click to show missing sample names.</a>
                <div id="goleft_empty_samples" style="display:none;">
                    <ul><li><code>{"</code></li>, <li><code>".join(bin_plot_data_empty_samples)}</code></li></ul>
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
            plot=scatter.plot(bin_plot_data, pconfig),
            content=extra,
        )
