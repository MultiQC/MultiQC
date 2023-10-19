""" MultiQC submodule to parse output from Picard QualityYieldMetrics """

import logging
from collections import OrderedDict

from multiqc import config

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
    for f in self.find_log_files(f"{self.anchor}/quality_yield_metrics", filehandles=True):
        s_name = f["s_name"]
        for line in f["f"]:
            maybe_s_name = self.extract_sample_name(line, f, picard_tool="CollectQualityYieldMetrics")
            if maybe_s_name:
                s_name = maybe_s_name

            if self.is_line_right_before_table(line, picard_class="QualityYieldMetrics"):
                if header != f["f"].readline().strip().split("\t"):
                    continue
                line = f["f"].readline()
                fields = [int(field) for field in line.strip("\n").split("\t")]
                if s_name in all_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                all_data[s_name] = dict(zip(header, fields))
                self.add_data_source(f, section="QualityYieldMetrics")

    # Filter to strip out ignored sample names
    all_data = self.ignore_samples(all_data)

    if not all_data:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Write parsed data to a file
    self.write_data_file(all_data, f"multiqc_{self.anchor}_QualityYieldMetrics")

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
