'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers data from _overall_mean_cov.csv files.
It relies on the following official source:
Overall Mean Coverage Report, page 196.
https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf
'''

import logging
import re
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule

from .utils import make_log_report

# Initialise the logger.
log = logging.getLogger(__name__)


class DragenOverallMeanCovMetrics(BaseMultiqcModule):
    """Public members of the DragenOverallMeanCovMetrics module are available to
    other dragen modules. Defines 2 public functions and a private data holder.
    """

    def collect_overall_mean_cov_data(self):
        """Collects raw data for coverage_metrics.py from overall_mean_cov.csv files."""
        overall_mean_cov_data = defaultdict(lambda: defaultdict(dict))

        for file in self.find_log_files("dragen/overall_mean_cov_metrics"):
            out = parse_overall_mean_cov(file)

            if out["success"]:
                self.add_data_source(file, section="overall_mean_cov_metrics")

                # Data for the coverage_metrics.py module.
                if out["phenotype"]:
                    sample, phenotype, data = out["sample"], out["phenotype"], out["data"]
                    overall_mean_cov_data[file["root"]][sample][phenotype] = data
                # Currently there is no need to support other files. Pass for now.
                else:
                    pass

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, file["s_name"])

        # Report found info/warnings/errors, which were collected while
        # calling the coverage_parser and constructing cov_headers.
        # You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "overall_mean_cov_metrics")

        # No need to write the data.
        self.write_data_file(overall_mean_cov_data, "dragen_overall_mean_cov_data")

        # Just to pass the pre-commit
        # doi=

        return overall_mean_cov_data


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-prefix>_overall_mean_cov<arbitrary-suffix>.csv

debug/info:
- unknown_metrics are those, which could not be recognized by the code.
  Such metrics will be ignored.
- unusual_values are those except for float.
'''
log_data = {
    "invalid_file_names": defaultdict(list),
    "unknown_metrics": [],
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
}

# Official structure of files:   _overall_mean_cov.csv
# Accepted structure of files: .+_overall_mean_cov.*.csv
GEN_FILE_RGX = re.compile(r"(.+)_overall_mean_cov(.*)\.csv$")

# Special case. Coverage metrics files have the following structure:
# <output-prefix>.<coverage-region-prefix>_overall_mean_cov<arbitrary-suffix>.csv
COV_FILE_RGX = re.compile(r"(.+)\.(.+)_overall_mean_cov(.*)\.csv$")

# General structure of lines is not defined.
# Currently only 1 metric is present in the standard. It substitutes the line's regex.
AVG_RGX = re.compile(r"Average alignment coverage over ([^,]+),([^,]+)$", re.IGNORECASE)


def parse_overall_mean_cov(file_handler):
    """Parser for _overall_mean_cov.csv files."""

    root = file_handler["root"]
    file = file_handler["fn"]

    # First check the general structure.
    file_match = GEN_FILE_RGX.search(file)
    if file_match:
        sample, phenotype = file_match.group(1), None
        file_match = COV_FILE_RGX.search(file)
        if file_match:
            sample, phenotype, suffix = file_match.groups()
            if suffix:
                # Only non-'tumor' strings are concatenated with the sample.
                if not re.search("tumor", suffix, re.IGNORECASE):
                    sample += suffix
    else:
        log_data["invalid_file_names"][root].append(file)
        return {"success": False}

    # There is still a possibility for failing(eg empty file, unknown metric).
    success = False
    source = None  # File name from "Average alignment coverage over file" metric.
    value = None

    # A loop is created in advance. Maybe new metrics will be added in the future.
    for line in file_handler["f"].splitlines():
        line_match = AVG_RGX.search(line)
        # If line is fine then extract the necessary fields.
        if line_match:
            success = True
            source, value = line_match.groups()

        # Otherwise check if line is empty. If not then report it and go to the next line.
        else:
            if not re.search(r"^\s*$", line):
                log_data["unknown_metrics"].append(line)
            continue

        # Try to convert the value. It shall be float.
        try:
            value = float(value)
        except ValueError:
            log_data["unusual_values"][root][file][line] = value

    return {
        "success": success,
        "sample": sample,
        "phenotype": phenotype,
        "data": {"source_file": source, "value": value},
    }
