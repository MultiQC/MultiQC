from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union

import logging

log = logging.getLogger(__name__)


def parse_reports(self) -> int:
    """Find and parse ngsbits SampleGender TSV output files."""

    samplegender_data: Dict[str, Dict[str, Union[float, int]]] = dict()
    for f in self.find_log_files("ngsbits/samplegender"):
        self.add_data_source(f)
        print(f["f"])  # File contents
        print(f["s_name"])  # Sample name (from cleaned filename)
        print(f["fn"])  # Filename
        print(f["root"])  # Directory file was in
        s_name = f["s_name"]
        if s_name in samplegender_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        samplegender_data[s_name] = parse_file(f["f"])
        print(samplegender_data)
        # Filter to strip out ignored sample names
        samplegender_data = self.ignore_samples(samplegender_data)

    n_reports_found = len(samplegender_data)
    if n_reports_found == 0:
        return 0

    # Write parsed report data to a file
    self.write_data_file(samplegender_data, "multiqc_ngsbits_samplegender")

    self.add_software_version(None)

    # Add SampleGender Table
    config_table = {
        "id": "ngsbits_samplegender",
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
        plot=table.plot(data=samplegender_data, headers=headers, pconfig=config_table),
    )

    self.general_stats_addcols(samplegender_data, headers)

    return len(samplegender_data)


def parse_file(f) -> Dict[str, Union[float, int, str]]:
    data = {}
    lines = f.splitlines()

    # Use the first line as the header (column names)
    headers = lines[0].strip().split("\t")

    # Process the second line, using the headers as key names
    values = lines[1].strip().split("\t")

    # Create a dictionary where headers[i] are the keys and values[i] are the corresponding values
    data = {
        headers[i]: float(values[i]) if i > 1 else values[i]  # Convert only numeric columns to float
        for i in range(1, len(values))
    }

    return data
