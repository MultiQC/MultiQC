from __future__ import absolute_import
from collections import OrderedDict, defaultdict

from multiqc.plots import linegraph
from .mapping_metrics import DragenMappingMetics
from .fragment_length import DragenFragmentLength
from .ploidy_estimation_metrics import DragenPloidyEstimationMetrics
from .time_metrics import DragenTimeMetrics
from .vc_metrics import DragenVCMetrics
from .coverage import DragenCoverage

import logging
log = logging.getLogger(__name__)


class MultiqcModule(DragenMappingMetics, DragenFragmentLength, DragenPloidyEstimationMetrics, DragenTimeMetrics, DragenVCMetrics, DragenCoverage):
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

        # collecting sample names based on mapping statistics file,
        # because the coverage stats like <tumor_name>.wgs_contig_mean_cov_normal.csv don't specify the normal sample name,
        self.sample_names_by_output_prefixes = defaultdict(dict)  # { <output-prefix>: {tumor: <tumor-name>, normal: <normal-name>} }

        self.parse_mapping_metrics()
        # T_SRR7890936_50pc.mapping_metrics.csv              - general stats table, barplots and beeswarm

        self.parse_fragment_length_hist()
        # T_SRR7890936_50pc.fragment_length_hist.csv         - histogram plot

        self.parse_ploidy_estimation_metrics()
        # T_SRR7890936_50pc.ploidy_estimation_metrics.csv    - add just PLOIDY ESTIMATION,,Ploidy estimation,X0  into gen stats

        self.parse_time_metrics()
        # T_SRR7890936_50pc.time_metrics.csv                 - perhaps a barplot?

        self.parse_vc_metrics()
        # T_SRR7890936_50pc.vc_metrics.csv                   - a table

        self.parse_coverage()
        # T_SRR7890936_50pc.wgs_contig_mean_cov_normal.csv   - histogram or a plot like in mosdepth, with each chrom in X axis
        # T_SRR7890936_50pc.wgs_contig_mean_cov_tumor.csv    - same
        # T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv  - metrics into gen stats table, own table, plus histogram for "PCT of genome with coverage"
        # T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv   - same
        # T_SRR7890936_50pc.wgs_fine_hist_normal.csv         - distribution plot
        # T_SRR7890936_50pc.wgs_fine_hist_tumor.csv          - same
        # T_SRR7890936_50pc.wgs_hist_normal.csv
        # T_SRR7890936_50pc.wgs_hist_tumor.csv
        # T_SRR7890936_50pc.wgs_overall_mean_cov_normal.csv  - more accurate cov value. replace the value from mapping_metrics
        # T_SRR7890936_50pc.wgs_overall_mean_cov_tumor.csv   - same
