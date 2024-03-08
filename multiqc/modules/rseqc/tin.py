""" MultiQC submodule to parse output from RSeQC tin.py
http://rseqc.sourceforge.net/#tin-py """

import csv
import logging


# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find RSeQC tin reports and parse their data"""

    self.tin_data = dict()

    for f in self.find_log_files("rseqc/tin", filehandles=True):
        # Parse contents
        try:
            reader = csv.DictReader(f["f"], delimiter="\t")
            contents = next(reader)
        except (csv.Error, StopIteration) as e:
            log.error(f"Could not parse file '{f['fn']}': {e}")
            continue

        s_name = self.clean_s_name(contents["Bam_file"], f)
        contents.pop("Bam_file")
        self.tin_data[s_name] = contents

        # Add file to multiqc_sources.txt
        self.add_data_source(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, s_name)

    # Filter to strip out ignored sample names
    self.tin_data = self.ignore_samples(self.tin_data)

    if len(self.tin_data) > 0:
        # Write to file
        self.write_data_file(self.tin_data, "multiqc_rseqc_tin")

        # Add to general stats table
        self.general_stats_headers["TIN(stdev)"] = {
            "title": "TIN stdev",
            "description": "Standard Deviation for the Transcript Integriry Number (TIN)",
            "max": 100,
            "min": 0,
            "scale": "Reds",
            "hidden": True,
        }
        self.general_stats_headers["TIN(median)"] = {
            "title": "TIN",
            "description": "Median Transcript Integriry Number (TIN), indicating the RNA integrity of a sample",
            "max": 100,
            "min": 0,
            "scale": "RdBu",
        }

        for s_name in self.tin_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name]["TIN(median)"] = self.tin_data[s_name]["TIN(median)"]
            self.general_stats_data[s_name]["TIN(stdev)"] = self.tin_data[s_name]["TIN(stdev)"]

    return len(self.tin_data)
