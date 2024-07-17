import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .bamPEFragmentSizeDistribution import bamPEFragmentSizeDistributionMixin

from .bamPEFragmentSizeTable import bamPEFragmentSizeTableMixin
from .estimateReadFiltering import EstimateReadFilteringMixin
from .plotCorrelation import plotCorrelationMixin
from .plotCoverage import plotCoverageMixin
from .plotEnrichment import PlotEnrichmentMixin
from .plotFingerprint import PlotFingerprintMixin
from .plotPCA import plotPCAMixin
from .plotProfile import plotProfileMixin

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
    """
    The module for deepTools parses a number of the text files that deepTools can produce. In particular, the following are supported:

    - `bamPEFragmentSize --table`
    - `bamPEFragmentSize --outRawFragmentLengths`
    - `estimateReadFiltering`
    - `plotCoverage ---outRawCounts` (as well as the content written normally to the console)
    - `plotEnrichment --outRawCounts`
    - `plotFingerprint --outQualityMetrics --outRawCounts`
    - `plotPCA --outFileNameData`
    - `plotCorrelation --outFileCorMatrix`
    - `plotProfile --outFileNameData`

    Please be aware that some tools (namely, `plotFingerprint --outRawCounts` and `plotCoverage --outRawCounts`) are only supported as of deepTools version 2.6. For earlier output from `plotCoverage --outRawCounts`, you can use `#'chr' 'start' 'end'` in `search_patterns.yaml` (see [here](http://multiqc.info/docs/#module-search-patterns) for more details). Also for these types of files, you may need to increase the maximum file size supported by MultiQC (`log_filesize_limit` in the MultiQC configuration file). You can find details regarding the configuration file location [here](http://multiqc.info/docs/#configuring-multiqc).

    Note that sample names are parsed from the text files themselves, they are not derived from file names.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="deepTools",
            anchor="deepTools",
            target="deepTools",
            href="http://deeptools.readthedocs.io",
            info="Tools to process and analyze deep sequencing data.",
            extra="deepTools addresses the challenge of handling the large amounts of data that are now routinely"
            "generated from DNA sequencing centers. deepTools contains useful modules to process the mapped "
            "reads data for multiple quality checks, creating **normalized coverage files** in standard bedGraph "
            "and bigWig file formats, that allow comparison between different files (for example, treatment and control). "
            "Finally, using such normalized and standardized files, deepTools can create many publication-ready "
            "**visualizations** to identify enrichments and for functional annotations of the genome.",
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
