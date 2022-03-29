#!/usr/bin/env python

""" MultiQC module to parse output from Picard """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Import the Picard submodules
from . import AlignmentSummaryMetrics
from . import BaseDistributionByCycleMetrics
from . import CollectIlluminaBasecallingMetrics
from . import CollectIlluminaLaneMetrics
from . import CrosscheckFingerprints
from . import ExtractIlluminaBarcodes
from . import GcBiasMetrics
from . import HsMetrics
from . import InsertSizeMetrics
from . import MarkDuplicates
from . import MarkIlluminaAdapters
from . import OxoGMetrics
from . import QualityByCycleMetrics
from . import QualityScoreDistributionMetrics
from . import QualityYieldMetrics
from . import RnaSeqMetrics
from . import RrbsSummaryMetrics
from . import TargetedPcrMetrics
from . import ValidateSamFile
from . import VariantCallingMetrics
from . import WgsMetrics


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Picard is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module output if
    logs are found."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Picard",
            anchor="picard",
            href="http://broadinstitute.github.io/picard/",
            info="is a set of Java command line tools for manipulating high-" "throughput sequencing data.",
            # No DOI to cite // doi=
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["AlignmentMetrics"] = AlignmentSummaryMetrics.parse_reports(self)
        if n["AlignmentMetrics"] > 0:
            log.info("Found {} AlignmentSummaryMetrics reports".format(n["AlignmentMetrics"]))

        n["BaseDistributionByCycleMetrics"] = BaseDistributionByCycleMetrics.parse_reports(self)
        if n["BaseDistributionByCycleMetrics"] > 0:
            log.info("Found {} BaseDistributionByCycleMetrics reports".format(n["BaseDistributionByCycleMetrics"]))

        n["CrosscheckFingerprints"] = CrosscheckFingerprints.parse_reports(self)
        if n["CrosscheckFingerprints"] > 0:
            log.info("Found {} CrosscheckFingerprints reports".format(n["CrosscheckFingerprints"]))

        n["GcBiasMetrics"] = GcBiasMetrics.parse_reports(self)
        if n["GcBiasMetrics"] > 0:
            log.info("Found {} GcBiasMetrics reports".format(n["GcBiasMetrics"]))

        n["HsMetrics"] = HsMetrics.parse_reports(self)
        if n["HsMetrics"] > 0:
            log.info("Found {} HsMetrics reports".format(n["HsMetrics"]))

        n["InsertSizeMetrics"] = InsertSizeMetrics.parse_reports(self)
        if n["InsertSizeMetrics"] > 0:
            log.info("Found {} InsertSizeMetrics reports".format(n["InsertSizeMetrics"]))

        n["MarkDuplicates"] = MarkDuplicates.parse_reports(self)
        if n["MarkDuplicates"] > 0:
            log.info("Found {} MarkDuplicates reports".format(n["MarkDuplicates"]))

        n["OxoGMetrics"] = OxoGMetrics.parse_reports(self)
        if n["OxoGMetrics"] > 0:
            log.info("Found {} OxoGMetrics reports".format(n["OxoGMetrics"]))

        n["QualityByCycleMetrics"] = QualityByCycleMetrics.parse_reports(self)
        if n["QualityByCycleMetrics"] > 0:
            log.info("Found {} QualityByCycleMetrics reports".format(n["QualityByCycleMetrics"]))

        n["QualityScoreDistributionMetrics"] = QualityScoreDistributionMetrics.parse_reports(self)
        if n["QualityScoreDistributionMetrics"] > 0:
            log.info("Found {} QualityScoreDistributionMetrics reports".format(n["QualityScoreDistributionMetrics"]))

        n["QualityYieldMetrics"] = QualityYieldMetrics.parse_reports(self)
        if n["QualityYieldMetrics"] > 0:
            log.info("Found {} QualityYieldMetrics reports".format(n["QualityYieldMetrics"]))

        n["RnaSeqMetrics"] = RnaSeqMetrics.parse_reports(self)
        if n["RnaSeqMetrics"] > 0:
            log.info("Found {} RnaSeqMetrics reports".format(n["RnaSeqMetrics"]))

        n["RrbsSummaryMetrics"] = RrbsSummaryMetrics.parse_reports(self)
        if n["RrbsSummaryMetrics"] > 0:
            log.info("Found {} RrbsSummaryMetrics reports".format(n["RrbsSummaryMetrics"]))

        n["TargetedPcrMetrics"] = TargetedPcrMetrics.parse_reports(self)
        if n["TargetedPcrMetrics"] > 0:
            log.info("Found {} TargetedPcrMetrics reports".format(n["TargetedPcrMetrics"]))

        n["VariantCallingMetrics"] = VariantCallingMetrics.parse_reports(self)
        if n["VariantCallingMetrics"] > 0:
            log.info("Found {} VariantCallingMetrics reports".format(n["VariantCallingMetrics"]))

        n["ValidateSamFile"] = ValidateSamFile.parse_reports(self)
        if n["ValidateSamFile"] > 0:
            log.info("Found {} ValidateSamFile reports".format(n["ValidateSamFile"]))

        n["WgsMetrics"] = WgsMetrics.parse_reports(self)
        if n["WgsMetrics"] > 0:
            log.info("Found {} WgsMetrics reports".format(n["WgsMetrics"]))

        n["CollectIlluminaBasecallingMetrics"] = CollectIlluminaBasecallingMetrics.parse_reports(self)
        if n["CollectIlluminaBasecallingMetrics"] > 0:
            log.info(
                "Found {} CollectIlluminaBasecallingMetrics reports".format(n["CollectIlluminaBasecallingMetrics"])
            )

        n["CollectIlluminaLaneMetrics"] = CollectIlluminaLaneMetrics.parse_reports(self)
        if n["CollectIlluminaLaneMetrics"] > 0:
            log.info("Found {} CollectIlluminaLaneMetrics reports".format(n["CollectIlluminaLaneMetrics"]))

        n["ExtractIlluminaBarcodes"] = ExtractIlluminaBarcodes.parse_reports(self)
        if n["ExtractIlluminaBarcodes"] > 0:
            log.info("Found {} ExtractIlluminaBarcodes reports".format(n["ExtractIlluminaBarcodes"]))

        n["MarkIlluminaAdapters"] = MarkIlluminaAdapters.parse_reports(self)
        if n["MarkIlluminaAdapters"] > 0:
            log.info("Found {} MarkIlluminaAdapters reports".format(n["MarkIlluminaAdapters"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    # Helper functions
    def multiply_hundred(self, val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
