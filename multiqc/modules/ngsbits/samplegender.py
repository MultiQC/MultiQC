from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict

import logging

log = logging.getLogger(__name__)


def parse_reports(self) -> int:
    """Find and parse ngsbits SampleGender TSV output files."""

    samplegender_data = dict()
    for f in self.find_log_files("ngsbits/samplegender"):
        print(f["f"])  # File contents
        print(f["s_name"])  # Sample name (from cleaned filename)
        print(f["fn"])  # Filename
        print(f["root"])  # Directory file was in
        parsed_data = samplegender_parse_reports(f["f"])
        if parsed_data:
            if f["s_name"] in samplegender_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f, section="samplegender")
            samplegender_data[f["s_name"]] = parsed_data
        print(parsed_data)
    # Filter to strip out ignored sample names
    #    samplegender_data = self.ignore_samples(samplegender_data)

    n_reports_found = len(samplegender_data)
    if n_reports_found == 0:
        return 0

    # Write parsed report data to a file
    self.write_data_file(samplegender_data, "multiqc_ngsbits_samplegender.txt")

    self.add_software_version(None)

    # Add SampleGender Table
    config_table = {
        "id": "ngsbits_samplegender",
        "namespace": "samplegender",
        "title": "ngs-bits: SampleGender",
    }

    headers = {
        "file": {
            "title": "File Name",
            "description": "The name of the BAM file analyzed.",
            "namespace": "ngsbits",
        },
        "gender": {
            "title": "Predicted Gender",
            "description": "The predicted gender based on chromosome read ratios.",
            "namespace": "ngsbits",
        },
        "reads_chry": {
            "title": "Reads on ChrY",
            "description": "The number of reads mapped to the Y chromosome.",
            "namespace": "ngsbits",
            "min": 0,
        },
        "reads_chrx": {
            "title": "Reads on ChrX",
            "description": "The number of reads mapped to the X chromosome.",
            "namespace": "ngsbits",
            "min": 0,
        },
        "ratio_chry_chrx": {
            "title": "ChrY/ChrX Ratio",
            "description": "The ratio of reads mapped to ChrY vs ChrX.",
            "namespace": "ngsbits",
            "min": 0,
            "format": "{:.4f}",
        },
    }

    self.add_section(
        name="ngs-bits SampleGender",
        anchor="samplegender",
        description="SampleGender determines the gender of a sample from the BAM/CRAM file.",
        plot=table.plot(samplegender_data, headers, config_table),
    )

    return len(samplegender_data)


def samplegender_parse_reports(f):
    data = dict()
    headers = None

    for line in f:
        if line.startswith("#"):  # Skip header or comment lines
            continue
        parts = line.strip().split("\t")
        if headers is None:
            headers = parts  # Use the first non-comment line as headers
            continue
        if len(parts) >= 5:  # Ensure there are at least five columns
            sample_data = {
                "file": parts[0],
                "gender": parts[1],
                "reads_chry": int(parts[2]),
                "reads_chrx": int(parts[3]),
                "ratio_chry_chrx": float(parts[4]),
            }
            data[parts[0]] = sample_data  # Use the file name as the key

    return data
