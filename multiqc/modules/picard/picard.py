import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Picard submodules, each one matching a picard tool
from . import (
    AlignmentSummaryMetrics,
    BaseDistributionByCycleMetrics,
    IlluminaBasecallingMetrics,
    IlluminaLaneMetrics,
    CrosscheckFingerprints,
    ExtractIlluminaBarcodes,
    GcBiasMetrics,
    HsMetrics,
    InsertSizeMetrics,
    MarkDuplicates,
    MarkIlluminaAdapters,
    OxoGMetrics,
    QualityByCycleMetrics,
    QualityScoreDistributionMetrics,
    QualityYieldMetrics,
    RnaSeqMetrics,
    RrbsSummaryMetrics,
    TargetedPcrMetrics,
    ValidateSamFile,
    VariantCallingMetrics,
    WgsMetrics,
)

TOOLS = [
    m.__name__.split(".")[-1]
    for m in (
        AlignmentSummaryMetrics,
        BaseDistributionByCycleMetrics,
        IlluminaBasecallingMetrics,
        IlluminaLaneMetrics,
        CrosscheckFingerprints,
        ExtractIlluminaBarcodes,
        GcBiasMetrics,
        HsMetrics,
        InsertSizeMetrics,
        MarkDuplicates,
        MarkIlluminaAdapters,
        OxoGMetrics,
        QualityByCycleMetrics,
        QualityScoreDistributionMetrics,
        QualityYieldMetrics,
        RnaSeqMetrics,
        RrbsSummaryMetrics,
        TargetedPcrMetrics,
        ValidateSamFile,
        VariantCallingMetrics,
        WgsMetrics,
    )
]

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands:

    - `AlignmentSummaryMetrics`
    - `BaseDistributionByCycle`
    - `CollectIlluminaBasecallingMetrics`
    - `CollectIlluminaLaneMetrics`
    - `CrosscheckFingerprints`
    - `ExtractIlluminaBarcodes`
    - `GcBiasMetrics`
    - `HsMetrics`
    - `InsertSizeMetrics`
    - `MarkDuplicates`
    - `MarkIlluminaAdapters`
    - `OxoGMetrics`
    - `QualityByCycleMetrics`
    - `QualityScoreDistributionMetrics`
    - `QualityYieldMetrics`
    - `RnaSeqMetrics`
    - `RrbsSummaryMetrics`
    - `ValidateSamFile`
    - `VariantCallingMetrics`
    - `WgsMetrics`

    #### Coverage Levels

    It's possible to customise the HsMetrics _"Target Bases 30X"_ coverage and
    WgsMetrics _"Fraction of Bases over 30X"_ that are
    shown in the general statistics table. This must correspond to field names in the
    picard report, such as `PCT_TARGET_BASES_2X` / `PCT_10X`. Any numbers not found in the
    reports will be ignored.

    The coverage levels available for HsMetrics are
    [typically](http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics)
    1, 2, 10, 20, 30, 40, 50 and 100X.

    The coverage levels available for WgsMetrics are
    [typically](http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics)
    1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90 and 100X.

    To customise this, add the following to your MultiQC config:

    ```yaml
    picard_config:
      general_stats_target_coverage:
        - 10
        - 50
    ```

    #### CrosscheckFingerprints

    In addition to adding a table of results, a `Crosschecks All Expected` column will be added to the General Statistics. If all comparisons for a sample were `Expected`, then the value of the field will be `True` and green. If not it will be `False` and Red.

    You can customize the columns show in the CrosscheckFingerprints table with the config keys `CrosscheckFingerprints_table_cols` and `CrosscheckFingerprints_table_cols_hidden`. For example:

    ```yaml
    picard_config:
      CrosscheckFingerprints_table_cols:
        - RESULT
        - LOD_SCORE
      CrosscheckFingerprints_table_cols_hidden:
        - LEFT_LANE
        - RIGHT_LANE
    ```

    The column names will be normalized, ex `LOD_SCORE -> Lod score`.

    Note that if `CALCULATE_TUMOR_AWARE_RESULTS` was set to true on the CLI for any of the CrosscheckFingerprints result files, then the `LOD_SCORE_TUMOR_NORMAL` and `LOD_SCORE_NORMAL_TUMOR` will be displayed.

    #### HsMetrics

    Note that the _Target Region Coverage_ plot is generated using the `PCT_TARGET_BASES_` table columns from the HsMetrics output (not immediately obvious when looking at the log files).

    You can customize the columns shown in the HsMetrics table with the config keys `HsMetrics_table_cols` and `HsMetrics_table_cols_hidden`. For example:

    ```yaml
    picard_config:
      HsMetrics_table_cols:
        - NEAR_BAIT_BASES
        - OFF_BAIT_BASES
        - ON_BAIT_BASES
      HsMetrics_table_cols_hidden:
        - MAX_TARGET_COVERAGE
        - MEAN_BAIT_COVERAGE
        - MEAN_TARGET_COVERAGE
    ```

    Only values listed in `HsMetrics_table_cols` will be included in the table.
    Anything listed in `HsMetrics_table_cols_hidden` will be hidden by default.

    A similar config is available for customising the HsMetrics columns in the General Stats table:

    ```yaml
    picard_config:
      HsMetrics_genstats_table_cols:
        - NEAR_BAIT_BASES
      HsMetrics_genstats_table_cols_hidden:
        - MAX_TARGET_COVERAGE
    ```

    #### InsertSizeMetrics

    By default, the insert size plot is smoothed to contain a maximum of 500 data points per sample.
    This is to prevent the MultiQC report from being very large with big datasets.
    If you would like to customise this value to get a better resolution you can set the following
    MultiQC config values, with the new maximum number of points:

    ```yaml
    picard_config:
      insertsize_smooth_points: 10000
    ```

    The plotted maximum insert size can be set with:

    ```yaml
    picard_config:
      insertsize_xmax: 10000
    ```

    #### MarkDuplicates

    If a `BAM` file contains multiple read groups, Picard MarkDuplicates generates a report
    with multiple metric lines, one for each "library".

    By default, MultiQC will sum the values for every library it finds and recompute the
    `PERCENT_DUPLICATION` and `ESTIMATED_LIBRARY_SIZE` fields, giving a single set of results
    for each `BAM` file.

    If instead you would prefer each library to be treated as a separate sample, you can do so
    by setting the following MultiQC config:

    ```yaml
    picard_config:
      markdups_merge_multiple_libraries: False
    ```

    This prevents the merge and recalculation and appends the library name to the sample name.

    This behaviour is present in MultiQC since version 1.9. Before this, only the metrics from the
    first library were taken and all others were ignored.

    #### ValidateSamFile Search Pattern

    Generally, Picard adds identifiable content to the output of function calls. This is not the case for ValidateSamFile. In order to identify logs the MultiQC Picard submodule `ValidateSamFile` will search for filenames that contain 'validatesamfile' or 'ValidateSamFile'. One can customise the used search pattern by overwriting the `picard/sam_file_validation` pattern in your MultiQC config. For example:

    ```yaml
    sp:
      picard/sam_file_validation:
        fn: "*[Vv]alidate[Ss]am[Ff]ile*"
    ```

    #### WgsMetrics

    The coverage histogram from Picard typically shows a normal distribution with a very long tail.
    To make the plot easier to view, by default the module plots the line up to 99% of the data.
    This typically removes the long tail and gives a more useful graph.

    If you would like, you can set a specific value for the maximum coverage to cut the graph at.
    By setting this to a very large value, you will disable the cutting (the graph will automatically
    limit the axis at the maximum data point). You can do this as follows:

    ```yaml
    picard_config:
      wgsmetrics_histogram_max_cov: 500
    ```

    If running with very high coverage samples or using the Picard `CAP_COVERAGE` option,
    the coverage histogram can become very large indeed. For eaxmple, if reporting coverages of 1 million,
    it will have 1 million data points per sample. That can crash the browser and take a long time to run.

    There are two customisation MultiQC options to help with this.
    Firstly, MultiQC will automatically "smooth" the histogram to a maximum of `1000` data points by binning.
    This should stop the browser from crashing. You can tweak how many bins are used with the following:

    ```yaml
    picard_config:
      wgsmetrics_histogram_smooth: 1000
    ```

    Change `1000` to whatever number you want. If you don't want any smoothing, set it to a very high number
    bigger than the number of data points you have.

    Secondly, if you would prefer to instead simply skip the histogram, you can set the following:

    ```yaml
    picard_config:
      wgsmetrics_skip_histogram: True
    ```

    This will omit that section from the report entirely, and also skip parsing the histogram data.
    By specifying this option you may speed up the run time for MultiQC with these types of files
    significantly.

    #### Sample names

    MultiQC supports outputs from multiple runs of a Picard tool merged together into one
    file. In order to handle multiple sample data in on file correctly, MultiQC needed
    to take the sample name elsewhere rather than the file name. For this reason, MultiQC
    attempts to parse the command line recorded in the output header. For example, an
    output from the `GcBias` tool contains a header line like this:

    ```
    # net.sf.picard.analysis.CollectGcBiasMetrics REFERENCE_SEQUENCE=/reference/genome.fa
    INPUT=/alignments/P0001_101/P0001_101.bam OUTPUT=P0001_101.collectGcBias.txt ...
    ```

    MultiQC would extract the BAM file name that goes after `INPUT=` and take `P0001_101`
    as a sample name. If MultiQC fails to parse the command line for any reason, it will
    fall back to using the file name. It is also possible to force using the file names
    as sample names by enabling the following config option:

    ```yaml
    picard_config:
      s_name_filenames: true
    ```
    """

    def __init__(
        self,
        name="Picard",
        anchor="picard",
        href="http://broadinstitute.github.io/picard/",
        info="Tools for manipulating high-throughput sequencing data.",
        # No DOI to cite // doi=
        tools=tuple(TOOLS),
    ):
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        self.samples_parsed_by_tool = dict()

        for tool in tools:
            log.debug(f"Running picard tool {tool}")
            mod = globals()[tool]
            func = getattr(mod, "parse_reports", None)
            if func is not None:
                self.samples_parsed_by_tool[tool] = func(self)
                if len(self.samples_parsed_by_tool[tool]) > 0:
                    log.info(f"Found {len(self.samples_parsed_by_tool[tool])} {tool} reports")

        # Exit if we didn't find anything
        if all(len(v) == 0 for v in self.samples_parsed_by_tool.values()):
            raise ModuleNoSamplesFound
