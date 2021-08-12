#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard QualityYieldMetrics """

from collections import OrderedDict, defaultdict
import logging
import os
import re

from multiqc import config
from multiqc.plots import table, linegraph
from .util import read_sample_name

# Initialise the logger
log = logging.getLogger(__name__)

DESC = OrderedDict(
    [
        ("TOTAL_READS", "The total number of reads in the input file"),
        ("PF_READS", "The number of reads that are PF - pass filter"),
        ("READ_LENGTH", "The average read length of all the reads (will be fixed for a lane)"),
        ("TOTAL_BASES", "The total number of bases in all reads"),
        ("PF_BASES", "The total number of bases in all PF reads"),
        ("Q20_BASES", "The number of bases in all reads that achieve quality score 20 or higher"),
        ("PF_Q20_BASES", "The number of bases in PF reads that achieve quality score 20 or higher"),
        ("Q30_BASES", "The number of bases in all reads that achieve quality score 30 or higher"),
        ("PF_Q30_BASES", "The number of bases in PF reads that achieve quality score 30 or higher"),
        ("Q20_EQUIVALENT_YIELD", "The sum of quality scores of all bases divided by 20"),
        ("PF_Q20_EQUIVALENT_YIELD", "The sum of quality scores of all bases in PF reads divided by 20"),
    ]
)


def parse_reports(self):
    """Find Picard QualityYieldMetrics reports and parse their data"""

    # Set up vars
    all_data = OrderedDict()
    header = list(DESC.keys())

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/quality_yield_metrics", filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        commadecimal = None

        lines = iter(f["f"])

        clean_fn = lambda n: self.clean_s_name(n, f)
        s_name = read_sample_name(lines, clean_fn, "CollectQualityYieldMetrics")

        if s_name is None:
            continue

        sample_data = dict()
        try:
            # skip to the histogram
            line = next(lines)
            while not line.startswith("## METRICS CLASS"):
                line = next(lines)

            # check the header
            line = next(lines)
            if header != line.strip().split("\t"):
                continue

            # one row
            line = next(lines)
            fields = [int(field) for field in line.strip("\n").split("\t")]
            sample_data = OrderedDict(zip(header, fields))
        except StopIteration:
            pass

        if sample_data:
            all_data[s_name] = OrderedDict(zip(header, fields))

    # Filter to strip out ignored sample names
    all_data = self.ignore_samples(all_data)

    if not all_data:
        return 0

    # Write parsed data to a file
    self.write_data_file(all_data, "multiqc_picard_QualityYieldMetrics")

    # Add to the general stats table
    headers = {
        "TOTAL_READS": {
            "title": "{} Reads".format(config.read_count_prefix),
            "description": "The total number of reads in the input file ({})".format(config.read_count_desc),
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Blues",
            "shared_key": "read_count",
        }
    }
    self.general_stats_addcols(all_data, headers)

    return len(all_data)
