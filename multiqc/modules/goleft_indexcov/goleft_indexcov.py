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
        cov_plot = self.coverage_plot()
        roc_plot = self.roc_plot()
        bin_plot = self.bin_plot()
        if not cov_plot:
            log.debug("Did not find goleft indexcov outputs in {}".format(config.analysis_dir))
            raise UserWarning

        self.sections = list()
        if cov_plot:
            self.sections.append({
                'name': 'Read coverage',
                'anchor': 'goleft_indexcov-coverage',
                'content': cov_plot
            })
        if roc_plot:
            self.sections.append({
                'name': 'Scaled coverage ROC plot',
                'anchor': 'goleft_indexcov-roc',
                'content': roc_plot
            })
        if bin_plot:
            self.sections.append({
                'name': 'Problem coverage bins',
                'anchor': 'goleft_indexcov-roc',
                'content': bin_plot
            })

    def _short_chrom(self, chrom):
        """Plot standard chromosomes + X, sorted numerically.
        """
        chrom = chrom.replace("chr", "")
        try:
            return int(chrom)
        except ValueError:
            if chrom in set(["X"]):
                return chrom

    def coverage_plot(self):
        max_y = 2.5
        data = collections.defaultdict(lambda: collections.defaultdict(dict))
        for fn in self.find_log_files(config.sp["goleft_indexcov"]["coverage"], filehandles=True):
            header = fn['f'].readline()
            sample_names = [self._clean_name(x) for x in header.strip().split()[3:]]
            for parts in (l.rstrip().split() for l in fn['f']):
                chrom, start, end = parts[:3]
                sample_covs = parts[3:]
                if self._short_chrom(chrom) is not None:
                    for cov, sample in zip(sample_covs, sample_names):
                        data[chrom][sample][int(end) + int(start) // 2] = min(max_y, float(cov))

        if data:
            chroms = sorted(data.keys(), key=self._short_chrom)
            pconfig = {
                'id': 'goleft_indexcov-coverage-plot',
                'title': 'Scaled coverage by chromosome',
                'ylab': 'Scaled coverage',
                'xlab': 'Position',
                'ymin': 0,
                'ymax': max_y,
                'xDecimals': False,
                'data_labels': [{"name": self._short_chrom(c)} for c in chroms]}
            return linegraph.plot([data[c] for c in chroms], pconfig)

    def roc_plot(self):
        data = collections.defaultdict(lambda: collections.defaultdict(dict))
        for fn in self.find_log_files(config.sp["goleft_indexcov"]["roc"], filehandles=True):
            header = fn['f'].readline()
            sample_names = [self._clean_name(x) for x in header.strip().split()[2:]]
            for parts in (l.rstrip().split() for l in fn['f']):
                chrom, cov = parts[:2]
                sample_vals = parts[2:]
                if self._short_chrom(chrom) is not None:
                    for val, sample in zip(sample_vals, sample_names):
                        data[chrom][sample][float(cov)] = float(val)
        if data:
            chroms = sorted(data.keys(), key=self._short_chrom)
            pconfig = {
                'id': 'goleft_indexcov-roc-plot',
                'title': 'ROC: genome coverage per scaled depth by chromosome',
                'xlab': 'Scaled coverage',
                'xlab': 'Proportion of regions covered',
                'ymin': 0, 'ymax': 1.0,
                'xmin': 0, 'xmax': 1.5,
                'data_labels': [{"name": self._short_chrom(c)} for c in chroms]}
            return linegraph.plot([data[c] for c in chroms], pconfig)

    def _clean_name(self, name):
        """Remove standard extra extensions from sample names.
        """
        return name.replace(".bam", "")

    def bin_plot(self):
        desc = '<p>This plot identifies problematic samples using binned coverage distributions. \n\
        We expect bins to be around 1, so deviations from this indicate problems. \n\
        Low coverage bins (< 0.15) on the x-axis have regions with low or missing coverage. \n\
        Higher values indicate truncated BAM files or missing data. \n\
        Bins with skewed distributions (<0.85 or >1.15) on the y-axis detect dosage bias. \n\
        Large values on the y-axis are likely to impact CNV and structural variant calling. \n\
        See the \n\
        <a href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md" target="_blank">goleft indexcov bin documentation</a> \n\
        for more details.</p>'

        data = {}
        for fn in self.find_log_files(config.sp["goleft_indexcov"]["ped"], filehandles=True):
            header = fn['f'].readline()[1:].strip().split("\t")
            for sample_parts in (l.split("\t") for l in fn['f']):
                cur = dict(zip(header, sample_parts))
                cur["sample_id"] = self._clean_name(cur["sample_id"])
                data[cur["sample_id"]] =  {"x": int(cur["bins.lo"]),
                                           "y": int(cur["bins.out"])}
        if data:
            pconfig = {
                'id': 'goleft_indexcov-bin-plot',
                'title': 'Problematic low and non-uniform coverage bins',
                'xDecimals': False, 'yDecimals': False,
                'xlab': 'Total bins with depth < 0.15',
                'ylab': 'Total bins with depth outside of (0.85, 1.15)'}
            return desc + scatter.plot(data, pconfig)
