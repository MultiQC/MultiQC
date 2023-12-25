"""MultiQC module to parse the output from deepTools"""
import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .bamPEFragmentSizeDistribution import bamPEFragmentSizeDistributionMixin

# deepTools modules
from .bamPEFragmentSizeTable import bamPEFragmentSizeTableMixin
from .estimateReadFiltering import EstimateReadFilteringMixin
from .plotCorrelation import plotCorrelationMixin
from .plotCoverage import plotCoverageMixin
from .plotEnrichment import PlotEnrichmentMixin
from .plotFingerprint import PlotFingerprintMixin
from .plotPCA import plotPCAMixin
from .plotProfile import plotProfileMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(
    BaseMultiqcModule,
    bamPEFragmentSizeTableMixin,
    bamPEFragmentSizeDistributionMixin,
    EstimateReadFilteringMixin,
    plotCoverageMixin,
    PlotEnrichmentMixin,
    PlotFingerprintMixin,
    plotProfileMixin,
    plotPCAMixin,
    plotCorrelationMixin,
):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="deepTools",
            anchor="deepTools",
            target="deepTools",
            href="http://deeptools.readthedocs.io",
            info=" is a suite of tools to process and analyze deep sequencing data.",
            doi="10.1093/nar/gkw257",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # plotCorrelation
        n["plotCorrelation"] = self.parse_plotCorrelation()
        if n["plotCorrelation"] > 0:
            log.debug(f"Found {n['plotCorrelation']} deepTools plotCorrelation samples")

        # plotPCA
        n["plotPCA"] = self.parse_plotPCA()
        if n["plotPCA"] > 0:
            log.debug(f"Found {n['plotPCA']} deepTools plotPCA samples")

        # plotEnrichment
        n["plotEnrichment"] = self.parse_plot_enrichment()
        if n["plotEnrichment"] > 0:
            log.debug(f"Found {n['plotEnrichment']} deepTools plotEnrichment samples")

        # plotFingerprint
        n["plotFingerprintOutQualityMetrics"], n["plotFingerprintOutRawCounts"] = self.parse_plotFingerprint()
        if n["plotFingerprintOutQualityMetrics"] + n["plotFingerprintOutRawCounts"] > 0:
            extra = ""
            if n["plotFingerprintOutRawCounts"] == 0:
                extra = (
                    " (you may need to increase the maximum log file size to find plotFingerprint --outRawCounts files)"
                )
            log.debug(
                "Found {} and {} deepTools plotFingerprint --outQualityMetrics and --outRawCounts samples, respectively{}".format(
                    n["plotFingerprintOutQualityMetrics"], n["plotFingerprintOutRawCounts"], extra
                )
            )

        # bamPEFragmentSizeDistribution
        n["bamPEFragmentSizeDistribution"] = self.parse_bamPEFragmentSizeDistribution()
        if n["bamPEFragmentSizeDistribution"] > 0:
            log.debug(
                "Found {} deepTools 'bamPEFragmentSize --outRawFragmentLengths' samples".format(
                    n["bamPEFragmentSizeDistribution"]
                )
            )

        # bamPEFragmentSizeTable
        n["bamPEFragmentSize"] = self.parse_bamPEFragmentSize()
        if n["bamPEFragmentSize"] > 0:
            log.debug(f"Found {n['bamPEFragmentSize']} deepTools 'bamPEFragmentSize --table' samples")

        # plotProfile
        n["plotProfile"] = self.parse_plotProfile()
        if n["plotProfile"] > 0:
            log.debug(f"Found {n['plotProfile']} deepTools plotProfile samples")

        # plotCoverage
        n["plotCoverageStdout"], n["plotCoverageOutRawCounts"] = self.parse_plotCoverage()
        if n["plotCoverageStdout"] + n["plotCoverageOutRawCounts"] > 0:
            extra = ""
            if n["plotCoverageOutRawCounts"] == 0:
                extra = (
                    " (you may need to increase the maximum log file size to find plotCoverage --outRawCounts files)"
                )
            log.debug(
                "Found {} and {} deepTools plotCoverage standard output and --outRawCounts samples, respectively{}".format(
                    n["plotCoverageStdout"], n["plotCoverageOutRawCounts"], extra
                )
            )

        # estimateReadFiltering
        n["estimateReadFiltering"] = self.parse_estimate_read_filtering()
        if n["estimateReadFiltering"] > 0:
            log.debug(f"Found {n['estimateReadFiltering']} deepTools estimateReadFiltering samples")

        tot = sum(n.values())
        if tot > 0:
            log.info(f"Found {tot} total deepTools samples")
        else:
            raise ModuleNoSamplesFound

    def _int(self, val):
        """Avoids Python3 error:
        ValueError: invalid literal for self._int() with base 10: '1.0'
        """
        return int(round(float(val)))
