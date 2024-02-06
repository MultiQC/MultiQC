""" MultiQC submodule to parse output from Picard QualityYieldMetrics """

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.picard import util

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


def parse_reports(module):
    """Find Picard QualityYieldMetrics reports and parse their data"""

    # Set up vars
    data_by_sample = dict()
    expected_header = list(DESC.keys())

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/quality_yield_metrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectQualityYieldMetrics",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="QualityYieldMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                if keys != expected_header:
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                module.add_data_source(f, s_name, section="QualityYieldMetrics")

                vals = []
                for v in f["f"].readline().strip("\n").split("\t"):
                    try:
                        v = int(v)
                    except ValueError:
                        pass
                    vals.append(v)
                data_by_sample[s_name] = dict(zip(keys, vals))
                s_name = None

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_QualityYieldMetrics")

    # Add to the general stats table
    headers = {
        "TOTAL_READS": {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"The total number of reads in the input file ({config.read_count_desc})",
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Blues",
            "shared_key": "read_count",
        }
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="QualityYieldMetrics")

    return len(data_by_sample)
