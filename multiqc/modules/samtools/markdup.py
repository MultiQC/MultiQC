# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from Samtools markdup """

import logging
import json
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MarkdupReportMixin:
    def parse_samtools_markdup(self):
        """Find Samtools markdup logs and parse their data"""

        self.samtools_markdup = dict()
        for f in self.find_log_files("samtools/markdup"):
            parsed_data = self.parse_report(f["f"])
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

        # Add to the general stats table
        general_stats_headers = {
            "PERCENT DUPLICATION": {
                "title": "Duplicates",
                "description": "markdup - Percent Duplication",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
            },
        }
        self.general_stats_addcols(self.samtools_markdup, general_stats_headers)

        # Prep the data for the plots
        # FEATURE PARITY with PICARD MarkDuplicates
        # UNIQUE PAIRS, UNIQUE UNPAIRED, DUPLICATE PAIR,DUPLICATE PAIR OPTICAL, DUPLICATE UNPAIRED, UNMAPPED (EXCLUDED)

        # Config for the plot
        pconfig = {
            "id": "samtools-markdup_deduplication",
            "title": "Samtools markdup: duplication Stats",
            "ylab": "# Reads",
        }

        self.add_section(
            name="Mark Duplicates",
            anchor="samtools-markdup",
            description="Number of reads, categorised by duplication state using Samtools.",
            plot=bargraph.plot(self.samtools_markdup, [], pconfig),
        )

        return len(self.samtools_markdup)

    def parse_report(self, f):
        """Parse the samtools markdup output file and return a dictionary of the data."""
        try:
            # check if the file is a json file
            parsed_data = json.loads(f)
        except ValueError:
            # if it is not a json file, then it is a text file
            parsed_data = {line.split(":")[0].strip(): line.split(":")[1].strip() for line in f.splitlines() if line}

        # Calculate the percent duplication
        parsed_data["PERCENT DUPLICATION"] = (
            int(parsed_data.get("DUPLICATE TOTAL", 0)) / int(parsed_data.get("READ")) * 100
        )

        # Double the "PAIR" data to fix the number of reads
        parsed_data["PAIRED"] = int(parsed_data.get("PAIRS", 0)) * 2
        parsed_data["DUPLICATE PAIR"] = int(parsed_data.get("DUPLICATE PAIR", 0)) * 2
        parsed_data["DUPLICATE PAIR OPTICAL"] = int(parsed_data.get("DUPLICATE PAIR OPTICAL", 0)) * 2

        # Drop unused data
        parsed_data.pop("COMMAND", None)
        parsed_data.pop("READ", None)
        parsed_data.pop("WRITTEN", None)
        parsed_data.pop("EXAMINED", None)
        parsed_data.pop("DUPLICATE TOTAL", None)
        parsed_data.pop("ESTIMATED_LIBRARY_SIZE", None)

        return parsed_data
