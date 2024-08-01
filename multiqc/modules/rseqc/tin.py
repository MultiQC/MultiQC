"""MultiQC submodule to parse output from RSeQC tin.py
http://rseqc.sourceforge.net/#tin-py"""

import csv
import logging
from typing import Dict

from multiqc import BaseMultiqcModule

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC tin reports and parse their data"""

    tin_data: Dict = dict()

    for f in module.find_log_files("rseqc/tin", filehandles=True):
        # Parse contents
        try:
            reader = csv.DictReader(f["f"], delimiter="\t")
            contents = next(reader)
        except (csv.Error, StopIteration) as e:
            log.error(f"Could not parse file '{f['fn']}': {e}")
            continue

        s_name = module.clean_s_name(contents["Bam_file"], f)
        contents.pop("Bam_file")
        tin_data[s_name] = contents

        # Add file to multiqc_sources.txt
        module.add_data_source(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        module.add_software_version(None, s_name)

    # Filter to strip out ignored sample names
    tin_data = module.ignore_samples(tin_data)

    if len(tin_data) > 0:
        # Write to file
        module.write_data_file(tin_data, "multiqc_rseqc_tin")

        # Add to general stats table
        headers = {
            "TIN(stdev)": {
                "title": "TIN stdev",
                "description": "Standard Deviation for the Transcript Integriry Number (TIN)",
                "max": 100,
                "min": 0,
                "scale": "Reds",
                "hidden": True,
            },
            "TIN(median)": {
                "title": "TIN",
                "description": "Median Transcript Integriry Number (TIN), indicating the RNA integrity of a sample",
                "max": 100,
                "min": 0,
                "scale": "RdBu",
            },
        }

        module.general_stats_addcols(tin_data, headers, namespace="TIN")

    return len(tin_data)
