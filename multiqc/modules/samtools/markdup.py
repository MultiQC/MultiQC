# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from Samtools markdup """

import logging
import json


# Initialise the logger
log = logging.getLogger(__name__)


class MarkdupReportMixin:
    def parse_samtools_markdup(self):
        """Find Samtools markdup logs and parse their data"""

        self.samtools_markdup = dict()
        for f in self.find_log_files("samtools/markdup"):
            parsed_data = parse_report(f["f"])
            print(parsed_data)
            if len(parsed_data) > 0:
                if f["s_name"] in self.samtools_markdup:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="markdup")
                self.samtools_markdup[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.samtools_markdup = self.ignore_samples(self.samtools_markdup)

        if len(self.samtools_markdup) == 0:
            return 0

        # Write parsed report data to a file (restructure first)
        self.write_data_file(self.samtools_markdup, "multiqc_samtools_markdup")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        return len(self.samtools_markdup)


def parse_report(f):
    """Parse the samtools markdup output file and return a dictionary of the data."""
    try:
        # check if the file is a json file
        parsed_data = json.loads(f)
    except ValueError:
        # if it is not a json file, then it is a text file
        parsed_data = {line.split(":")[0].strip(): line.split(":")[1].strip() for line in f.splitlines() if line}
    return parsed_data
