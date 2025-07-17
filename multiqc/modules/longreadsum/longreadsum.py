import copy
import logging
import os
import re
import json

import yaml

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, table
from multiqc.plots.table_object import ColumnMeta

from typing import Dict, Union

log = logging.getLogger(__name__)

# https://docs.seqera.io/multiqc/reports
# https://github.com/MultiQC/MultiQC/blob/main/.github/CONTRIBUTING.md
# https://docs.seqera.io/multiqc/development/modules


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses data in the `summary.txt` LongReadSum output files.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="LongReadSum",
            anchor="longreadsum",
            href="https://github.com/WGLab/LongReadSum.git",
            info="Quality control for long-read sequencing data.",
            extra="""
             LongReadSum: A fast and flexible quality control and signal
             summarization tool for long-read sequencing data.
            """,
            doi="10.1016/j.csbj.2025.01.019",
        )

        # Get data by sample
        data_by_sample: Dict[str, Dict[str, Union[float, int]]] = dict()
        for f in self.find_log_files("longreadsum/summary"):
            sample_name = f["s_name"]
            # log.debug(f"Processing file: {f['f']} for sample: {sample_name}")
            log.debug("Processing file for sample: %s", sample_name)
            if sample_name in data_by_sample:
                log.debug("Duplicate sample name found! Overwriting: %s", sample_name)

            data_by_sample[sample_name] = parse_file(f["f"], sample_name)

            # Add version information
            if "longreadsum_version" in data_by_sample[sample_name]:
                version = data_by_sample[sample_name]["longreadsum_version"]
                self.add_software_version(version, sample=sample_name)
                log.debug("Found LongReadSum version %s for sample %s", version, sample_name)
                # Remove the version from the data, not needed in the table
                del data_by_sample[sample_name]["longreadsum_version"]

        # Remove samples to ignore
        data_by_sample = self.ignore_samples(data_by_sample)

        # If no data found, raise an error
        if len(data_by_sample) == 0:
            log.warning("No LongReadSum summary files found.")
            raise ModuleNoSamplesFound(
                "No LongReadSum summary files found. "
                "Please ensure the output directory is correct and contains the expected files."
            )

        # Basic statistics table
        log.debug("Creating basic statistics table.")
        # self.general_stats_table(data_by_sample)
        self.general_stats_addcols(data_by_sample)


def parse_file(raw_text: str, sample: str) -> Dict[str, Union[float, int]]:
    """
    Parse the summary.json content from a raw text string and return a dictionary of values.
    """
    data = {}
    json_data = json.loads(raw_text)
    for key, value in json_data.items():
        data[key] = value

    return data
