from __future__ import absolute_import
from collections import OrderedDict, defaultdict

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


class MultiqcModule(DragenMappingMetics, DragenFragmentLength, DragenPloidyEstimationMetrics, DragenTimeMetrics, DragenVCMetrics,
                    DragenCoveragePerContig, DragenCoverageMetrics, DragenCoverageHist):
    """ Dragen has a number of different pipelines and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Dragen', anchor='Dragen', target='Dragen',
            href='https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html',
            info=(" is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data"
                  " using field-programmable gate array technology (FPGA)."))

        self.parse_coverage_hist()
        # <output prefix>.wgs_fine_hist_normal.csv         - distribution plot
        # <output prefix>.wgs_fine_hist_tumor.csv          - same

        self.parse_mapping_metrics()
        # <output prefix>.mapping_metrics.csv              - general stats table, barplots and beeswarm

        self.parse_coverage_metrics()
        # <output prefix>.wgs_coverage_metrics_normal.csv  - metrics into gen stats table, own table, plus histogram for "PCT of genome with coverage"
        # <output prefix>.wgs_coverage_metrics_tumor.csv   - same

        self.parse_coverage_per_contig()
        # <output prefix>.wgs_contig_mean_cov_normal.csv   - histogram or a plot like in mosdepth, with each chrom in X axis
        # <output prefix>.wgs_contig_mean_cov_tumor.csv    - same

        self.parse_fragment_length_hist()
        # <output prefix>.fragment_length_hist.csv         - histogram plot

        self.parse_ploidy_estimation_metrics()
        # <output prefix>.ploidy_estimation_metrics.csv    - add just PLOIDY ESTIMATION,,Ploidy estimation,X0  into gen stats

        self.parse_vc_metrics()
        # <output prefix>.vc_metrics.csv                   - a table

        self.parse_time_metrics()
        # <output prefix>.time_metrics.csv                 - perhaps a barplot?

