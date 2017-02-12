"""MultiQC module to plot output from goleft indexcov

https://github.com/brentp/goleft/tree/master/indexcov
"""
from __future__ import print_function
import collections
import logging
import os

from multiqc import config
from multiqc.plots import linegraph
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
        if not cov_plot and not roc_plot:
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
            sample = header.strip().split()[-1]
            for chrom, start, end, cov in (l.split() for l in fn['f']):
                if self._short_chrom(chrom) is not None:
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
            sample = header.strip().split()[-1]
            for chrom, cov, val in (l.split() for l in fn['f']):
                if self._short_chrom(chrom) is not None:
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
