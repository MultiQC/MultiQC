import logging

from ..base_module import ModuleNoSamplesFound
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
    """DRAGEN provides a number of different pipelines and outputs, including base calling, DNA and RNA alignment,
    post-alignment processing and variant calling, covering virtually all stages of typical NGS data processing.
    However, it can be treated as a fast aligner with additional features on top, as users will unlikely use any
    features without enabling DRAGEN mapping. So we will treat this module as an alignment tool module and
    place it accordingly in the module_order list, in docs, etc.

    The QC metrics DRAGEN generates resemble those of samtools-stats, qualimap, mosdepth, bcftools-stats and alike.
    Whenever possible, the visual output is made similar to those modules.

    Note that this MultiQC module supports some of DRAGEN output but not all. Contributions are welcome!

    The code is structured in a way so every mix-in parses one type of QC file that DRAGEN generates
    (e.g. *.mapping_metrics.csv, *.wgs_fine_hist_normal.csv, etc.). If a corresponding file is found, a mix-in adds
    a section into the report.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="DRAGEN",
            anchor="DRAGEN",
            target="DRAGEN",
            href="https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html",
            info=""" is a Bio-IT Platform that provides ultra-rapid secondary analysis of sequencing data
                     using field-programmable gate array technology (FPGA).""",
            # Can't find a DOI // doi=
        )

        samples_found = set()
        samples_found |= self.add_mapping_metrics()
        # <output prefix>.mapping_metrics.csv              - general stats table, a dedicated table, and a few barplots

        samples_found |= self.add_vc_metrics()
        # <output prefix>.vc_metrics.csv                   - a dedicated table and the total number of Variants into the general stats table

        samples_found |= self.add_ploidy_estimation_metrics()
        # <output prefix>.ploidy_estimation_metrics.csv    - add just Ploidy estimation into gen stats

        overall_mean_cov_data = self.collect_overall_mean_cov_data()
        # <output prefix>.<coverage region prefix>_overall_mean_cov<arbitrary suffix>.csv
        # This data will be used by in the DragenCoverageMetrics.

        samples_found |= self.add_coverage_metrics(overall_mean_cov_data)
        # <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv

        samples_found |= self.add_coverage_hist()

        # <output prefix>.wgs_fine_hist_normal.csv         - coverage distribution and cumulative coverage plots
        # <output prefix>.wgs_fine_hist_tumor.csv          - same

        samples_found |= self.add_coverage_per_contig()
        # <output prefix>.wgs_contig_mean_cov_normal.csv   - a histogram like in mosdepth, with each chrom as a category on X axis, plus a category for autosomal chromosomes average
        # <output prefix>.wgs_contig_mean_cov_tumor.csv    - same

        samples_found |= self.add_fragment_length_hist()
        # <output prefix>.fragment_length_hist.csv         - a histogram plot

        samples_found |= self.add_gc_metrics_hist()
        # <output prefix>.gc_metrics.csv

        samples_found |= self.add_trimmer_metrics()
        # <output prefix>.trimmer_metrics.csv

        samples_found |= self.add_time_metrics()
        # <output prefix>.time_metrics.csv

        samples_found |= self.add_rna_metrics()
        # <output prefix>.quant.metrics.csv

        samples_found |= self.add_rna_transcript_coverage()
        # <output prefix>.quant.transcript_coverage.txt

        samples_found |= self.add_sc_rna_metrics()
        # <output prefix>.scRNA.metrics.csv or <output prefix>.scRNA_metrics.csv

        samples_found |= self.add_sc_atac_metrics()
        # <output prefix>.scATAC.metrics.csv or <output prefix>.scATAC_metrics.csv

        if len(samples_found) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found samples: {len(samples_found)}")
