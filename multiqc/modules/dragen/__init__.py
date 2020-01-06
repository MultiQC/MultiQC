from __future__ import absolute_import

from .mapping_metrics import DragenMappingMetics
from .fragment_length import DragenFragmentLength
from .ploidy_estimation_metrics import DragenPloidyEstimationMetrics
from .time_metrics import DragenTimeMetrics
from .vc_metrics import DragenVCMetrics
from .coverage_per_contig import DragenCoveragePerContig
from .coverage_metrics import DragenCoverageMetrics
from .coverage_hist import DragenCoverageHist

import logging
log = logging.getLogger(__name__)


class MultiqcModule(DragenMappingMetics, DragenFragmentLength, DragenPloidyEstimationMetrics, DragenTimeMetrics,
                    DragenVCMetrics, DragenCoveragePerContig, DragenCoverageMetrics, DragenCoverageHist):
    """ Dragen has a number of differrent pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    For each stage, it generates QC files with metrics resembling those of samtools-stats, mosdepth, bcftools-stats
    and alike. This MultiQC module supports some of the output but not all.

    The code for each type of QC output is split into its own file, where a class mix-in is defined, that is
    inherited by this main MultiqcModule class. Each mix-in adds a section to the module output if corresponding
    QC files are found. """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Dragen', anchor='Dragen', target='Dragen',
            href='https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html',
            info=(" is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data"
                  " using field-programmable gate array technology (FPGA)."))

        self.parse_coverage_hist()
        # <output prefix>.wgs_fine_hist_normal.csv         - coverage distribution and cumulative coverage plots
        # <output prefix>.wgs_fine_hist_tumor.csv          - same

        self.parse_mapping_metrics()
        # <output prefix>.mapping_metrics.csv              - general stats table, a dedicated table, and a few barplots

        self.parse_coverage_metrics()
        # <output prefix>.wgs_coverage_metrics_normal.csv  - general stats table and a dedicated table
        # <output prefix>.wgs_coverage_metrics_tumor.csv   - same

        self.parse_coverage_per_contig()
        # <output prefix>.wgs_contig_mean_cov_normal.csv   - a histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average
        # <output prefix>.wgs_contig_mean_cov_tumor.csv    - same

        self.parse_fragment_length_hist()
        # <output prefix>.fragment_length_hist.csv         - a histogram plot

        self.parse_ploidy_estimation_metrics()
        # <output prefix>.ploidy_estimation_metrics.csv    - add just Ploidy estimation into gen stats

        self.parse_vc_metrics()
        # <output prefix>.vc_metrics.csv                   - a dedicated table and the total number of Variants into the general stats table

        self.parse_time_metrics()
        # <output prefix>.time_metrics.csv                 - ideally an overlapping barplot, but for now a beeswarm plot with dots for time each stage for each sample has finished

