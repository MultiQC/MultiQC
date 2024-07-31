import logging
from collections import defaultdict
from typing import Dict

from multiqc.base_module import ModuleNoSamplesFound
from .coverage_hist import DragenCoverageHist
from .coverage_metrics import DragenCoverageMetrics
from .coverage_per_contig import DragenCoveragePerContig
from .dragen_gc_metrics import DragenGcMetrics
from .fragment_length import DragenFragmentLength
from .mapping_metrics import DragenMappingMetics
from .overall_mean_cov import DragenOverallMeanCovMetrics
from .ploidy_estimation_metrics import DragenPloidyEstimationMetrics
from .rna_quant_metrics import DragenRnaQuantMetrics
from .rna_transcript_cov import DragenRnaTranscriptCoverage
from .sc_atac_metrics import DragenScAtacMetrics
from .sc_rna_metrics import DragenScRnaMetrics
from .time_metrics import DragenTimeMetrics
from .trimmer_metrics import DragenTrimmerMetrics
from .vc_metrics import DragenVCMetrics

log = logging.getLogger(__name__)


class MultiqcModule(
    DragenMappingMetics,
    DragenFragmentLength,
    DragenPloidyEstimationMetrics,
    DragenVCMetrics,
    DragenCoveragePerContig,
    DragenOverallMeanCovMetrics,
    DragenCoverageMetrics,
    DragenCoverageHist,
    DragenGcMetrics,
    DragenTrimmerMetrics,
    DragenTimeMetrics,
    DragenRnaQuantMetrics,
    DragenRnaTranscriptCoverage,
    DragenScRnaMetrics,
    DragenScAtacMetrics,
):
    """
    DRAGEN has a number of different pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    For each stage, it generates QC files with metrics resembling those of samtools-stats, mosdepth, bcftools-stats
    and alike. This MultiQC module supports some of the output but not all. Contributions are welcome!

    - `<output prefix>.wgs_fine_hist_<tumor|normal>.csv`
      - Coverage distribution and cumulative coverage plots
    - `<output prefix>.mapping_metrics.csv`
      - General stats table, a dedicated table, and a few barplots
    - `<output prefix>.wgs_coverage_metrics_<tumor|normal>.csv`
      - General stats table and a dedicated table
    - `<output prefix>.qc-coverage-region-<1|2|3>_coverage_metrics.csv`
      - General stats table and a dedicated table
    - `<output prefix>.wgs_contig_mean_cov_<tumor|normal>.csv`
      - A histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average
    - `<output prefix>.fragment_length_hist.csv`
      - A histogram plot
    - `<output prefix>.ploidy_estimation_metrics.csv`
      - Add just Ploidy estimation into the general stats table
    - `<output prefix>.vc_metrics.csv`
      - A dedicated table and the total number of Variants into the general stats table
    - `<output prefix>.gc_metrics.csv`
      - A histogram and summary statistics table on GC content metrics
    - `<output prefix>.trimmer_metrics.csv`
      - A summary table of tirmmer metrics
    - `<output prefix>.time_metrics.metrics`
      - A bar graph of the total run time and a breakdown of the run time of each individual step
    - `<output prefix>.quant.metrics.csv`
      - A bar graph of RNA fragments
    - `<output prefix>.quant.transcript_coverage.txt`
      - A line plot of average coverage along RNA transcripts
    - `<output prefix>.scRNA.metrics.csv` or `<output prefix>.scRNA_metrics.csv`
      - Summary table for single-cell RNA metrics
    - `<output prefix>.scATAC.metrics.csv` or `<output prefix>.scATAC_metrics.csv`
      - Summary table for single-cell ATAC metrics

    The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
    (e.g. *.mapping_metrics.csv, *.wgs_fine_hist_normal.csv, etc.). If a corresponding file is found, a mix-in adds
    a section into the report.

    DRAGEN can be treated as a fast aligner with additional features on top, as users will unlikely use any
    features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
    place it accordingly in the module_order list, in docs, etc.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="DRAGEN",
            anchor="DRAGEN",
            href="https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html",
            info="Illumina Bio-IT Platform that uses FPGA for secondary analysis of sequencing data.",
            # Can't find a DOI // doi=
        )

        functions = [
            self.add_mapping_metrics,
            # <output prefix>.mapping_metrics.csv           - general stats table, a dedicated table, and a few barplots
            self.add_vc_metrics,
            # <output prefix>.vc_metrics.csv                - a dedicated table and the total number of Variants into the general stats table
            self.add_ploidy_estimation_metrics,
            # <output prefix>.ploidy_estimation_metrics.csv - add just "Ploidy estimation" into gen stats
            self.collect_overall_mean_cov_data,
            # <output prefix>.<coverage region prefix>_overall_mean_cov<arbitrary suffix>.csv
            # This data will be also used by the DragenCoverageMetrics.add_coverage_metrics
            self.add_coverage_metrics,
            # <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv
            self.add_coverage_hist,
            # <output prefix>.wgs_fine_hist_normal.csv         - coverage distribution and cumulative coverage plots
            # <output prefix>.wgs_fine_hist_tumor.csv          - same
            self.add_coverage_per_contig,
            # <output prefix>.wgs_contig_mean_cov_normal.csv   - a histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average
            # <output prefix>.wgs_contig_mean_cov_tumor.csv    - same
            self.add_fragment_length_hist,
            # <output prefix>.fragment_length_hist.csv         - a histogram plot
            self.add_gc_metrics_hist,
            # <output prefix>.gc_metrics.csv
            self.add_trimmer_metrics,
            # <output prefix>.trimmer_metrics.csv
            self.add_time_metrics,
            # <output prefix>.time_metrics.csv
            self.add_rna_metrics,
            # <output prefix>.quant.metrics.csv
            self.add_rna_transcript_coverage,
            # <output prefix>.quant.transcript_coverage.txt
            self.add_sc_rna_metrics,
            # <output prefix>.scRNA.metrics.csv or <output prefix>.scRNA_metrics.csv
            self.add_sc_atac_metrics,
            # <output prefix>.scATAC.metrics.csv or <output prefix>.scATAC_metrics.csv
        ]

        # Populated by overall_mean_cov_data and used by add_coverage_hist
        self.overall_mean_cov_data: Dict[str, Dict[str, Dict]] = defaultdict(lambda: defaultdict(dict))

        self.samples_parsed_by_tool = dict()
        for func in functions:
            tool = func.__name__
            self.samples_parsed_by_tool[tool] = func()
            log.info(f"Found {len(self.samples_parsed_by_tool[tool])} {tool} reports")

        # Exit if we didn't find anything
        if all(len(v) == 0 for v in self.samples_parsed_by_tool.values()):
            raise ModuleNoSamplesFound
