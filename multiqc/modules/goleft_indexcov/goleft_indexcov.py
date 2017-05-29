"""MultiQC module to plot output from goleft indexcov

https://github.com/brentp/goleft/tree/master/indexcov
"""
from __future__ import print_function
import collections
import logging

from multiqc import config
from multiqc.plots import linegraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name='goleft indexcov', anchor='goleft_indexcov',
                                            href='https://github.com/brentp/goleft/tree/master/indexcov',
                                            info="quickly estimates coverage from a whole-genome bam index.")
        roc_plot = self.roc_plot()
        bin_plot = self.bin_plot()
        if not roc_plot and not bin_plot:
            log.debug("Did not find goleft indexcov outputs in {}".format(config.analysis_dir))
            raise UserWarning

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

    def roc_plot(self):
        desc = 'Coverage (ROC) plot that shows genome coverage at at given (scaled) depth. \n\
        Lower coverage samples have shorter curves where the proportion of regions covered \n\
        drops off more quickly. This indicates a higher fraction of low coverage regions.'
        max_chroms = 50
        data = collections.defaultdict(lambda: collections.defaultdict(dict))
        for fn in self.find_log_files('goleft_indexcov/roc', filehandles=True):
            header = fn['f'].readline()
            sample_names = [self.clean_s_name(x, fn["root"]) for x in header.strip().split()[2:]]
            for parts in (l.rstrip().split() for l in fn['f']):
                if len(parts) > 2:
                    chrom, cov = parts[:2]
                    sample_vals = parts[2:]
                    if self._short_chrom(chrom) is not None:
                        for val, sample in zip(sample_vals, sample_names):
                            data[chrom][sample][float(cov)] = float(val)

        # Filter to strip out ignored sample names
        for chrom in data:
            data[chrom] = self.ignore_samples(data[chrom])

        if data:
            def to_padded_str(x):
                x = self._short_chrom(x)
                try:
                    return "%06d" % x
                except TypeError:
                    return x
            chroms = sorted(data.keys(), key=to_padded_str)
            log.info("Found goleft indexcov ROC reports for %s samples" % (len(data[chroms[0]])))
            if len(chroms) > max_chroms:
                log.info("Too many chromosomes found: %s, limiting to %s" % (len(chroms), max_chroms))
                chroms = chroms[:max_chroms]
            pconfig = {
                'id': 'goleft_indexcov-roc-plot',
                'title': 'ROC: genome coverage per scaled depth by chromosome',
                'xlab': 'Scaled coverage',
                'ylab': 'Proportion of regions covered',
                'ymin': 0, 'ymax': 1.0,
                'xmin': 0, 'xmax': 1.5,
                'data_labels': [{"name": self._short_chrom(c)} for c in chroms]}
            self.add_section (
                name = 'Scaled coverage ROC plot',
                anchor = 'goleft_indexcov-roc',
                description = desc,
                plot = linegraph.plot([data[c] for c in chroms], pconfig)
            )
            return True
        else:
            return False

    def bin_plot(self):
        desc = 'This plot identifies problematic samples using binned coverage distributions. \n\
        We expect bins to be around 1, so deviations from this indicate problems. \n\
        Low coverage bins (< 0.15) on the x-axis have regions with low or missing coverage. \n\
        Higher values indicate truncated BAM files or missing data. \n\
        Bins with skewed distributions (<0.85 or >1.15) on the y-axis detect dosage bias. \n\
        Large values on the y-axis are likely to impact CNV and structural variant calling. \n\
        See the \n\
        <a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md" target="_blank">goleft indexcov bin documentation</a> \n\
        for more details.'

        data = {}
        for fn in self.find_log_files('goleft_indexcov/ped', filehandles=True):
            header = fn['f'].readline()[1:].strip().split("\t")
            for sample_parts in (l.split("\t") for l in fn['f']):
                cur = dict(zip(header, sample_parts))
                cur["sample_id"] = self.clean_s_name(cur["sample_id"], fn["root"])
                total = float(cur["bins.in"]) + float(cur["bins.out"])
                data[cur["sample_id"]] = {"x": float(cur["bins.lo"]) / total,
                                          "y": float(cur["bins.out"]) / total}

        # Filter to strip out ignored sample names
        data = self.ignore_samples(data)

        if data:
            log.info("Found goleft indexcov bin reports for %s samples" % (len(data)))
            pconfig = {
                'id': 'goleft_indexcov-bin-plot',
                'title': 'Problematic low and non-uniform coverage bins',
                'xlab': 'Proportion of bins with depth < 0.15',
                'ylab': 'Proportion of bins with depth outside of (0.85, 1.15)',
                'yCeiling': 1.0, 'yFloor': 0.0, 'xCeiling': 1.0, 'xFloor': 0.0}
            self.add_section (
                name = 'Problem coverage bins',
                anchor = 'goleft_indexcov-roc',
                description = desc,
                plot = scatter.plot(data, pconfig)
            )
            return True
        else:
            return False
