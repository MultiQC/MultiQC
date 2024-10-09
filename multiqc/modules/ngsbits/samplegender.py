from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union

import logging

log = logging.getLogger(__name__)


def parse_reports(self: BaseMultiqcModule) -> int:
    """Find and parse ngsbits SampleGender TSV output files."""

    samplegender_data: Dict[str, Dict[str, Union[float, str]]] = dict()
    for f in self.find_log_files("ngsbits/samplegender"):
        self.add_data_source(f)
        s_name = f["s_name"]
        s_name = s_name.replace("_ngsbits_sex", "")

        if s_name in samplegender_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        samplegender_data[s_name] = parse_file(f["f"])
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
        "id": "samplegender",
        "title": "SampleGender",
    }

    headers = {
        "gender": {
            "title": "Predicted Gender",
            "description": "The predicted gender based on chromosome read ratios.",
            "namespace": "ngsbits",
            "scale": False,
        },
        "reads_chry": {
            "title": "Reads on ChrY",
            "description": "The number of reads mapped to the Y chromosome.",
            "namespace": "ngsbits",
            "min": 0,
            "scale": "Blues",
            "format": "{:,.0f}",
        },
        "reads_chrx": {
            "title": "Reads on ChrX",
            "description": "The number of reads mapped to the X chromosome.",
            "namespace": "ngsbits",
            "min": 0,
            "scale": "Reds",
            "format": "{:,.0f}",
        },
        "ratio_chry_chrx": {
            "title": "ChrY/ChrX Ratio",
            "description": "The ratio of reads mapped to ChrY vs ChrX.",
            "namespace": "ngsbits",
            "min": 0,
            "format": "{:.4f}",
            "scale": "Purples",
        },
    }

    self.add_section(
        name="SampleGender",
        anchor="ngsbits-samplegender",
        description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/SampleGender.md" target="_blank">SampleGender</a>'
        " determines the gender of a sample from the BAM/CRAM file.",
        plot=table.plot(data=samplegender_data, headers=headers, pconfig=config_table),
    )

    for header in headers.values():
        header["hidden"] = True
    headers["gender"]["hidden"] = False

    self.general_stats_addcols(samplegender_data, headers)

    return len(samplegender_data)


def parse_file(f: str) -> Dict[str, Union[float, str]]:
    """
    Parses a single SampleGender TSV file content and returns a dictionary
    with the relevant data from columns 2-5.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()

    if len(lines) < 2:
        # Not enough data, return an empty dictionary
        return parsed_data

    headers = lines[0].strip().split("\t")[1:5]
    values = lines[1].strip().split("\t")[1:5]

    for key, value in zip(headers, values):
        try:
            parsed_data[key] = float(value)
        except ValueError:
            parsed_data[key] = value

    return parsed_data
