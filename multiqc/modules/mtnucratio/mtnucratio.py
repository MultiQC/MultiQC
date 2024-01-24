""" MultiQC module to parse output from mtnucratio """


import json
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """mtnucratio module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="mtnucratio",
            anchor="mtnucratio",
            href="http://www.github.com/apeltzer/MTNucRatioCalculator",
            info="is a tool to compute mt/nuc ratios for NGS datasets.",
            doi="10.1186/s13059-016-0918-z",
        )

        # Find and load any MTNUCRATIO reports
        self.mtnuc_data = dict()

        for f in self.find_log_files("mtnucratio", filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.mtnuc_data = self.ignore_samples(self.mtnuc_data)

        if len(self.mtnuc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.mtnuc_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.mtnuc_data, "multiqc_mtnucratio")

        # Basic Stats Table
        self.mtnucratio_general_stats_table()

    # Parse our nice little JSON file
    def parseJSON(self, f):
        """Parse the JSON output from mtnucratio and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
            if "metrics" not in parsed_json and "metadata" not in parsed_json:
                log.warning(f"No MTNUCRATIO JSON: '{f['fn']}'")
                return None
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse mtnucratio JSON: '{f['fn']}'")
            log.debug(e)
            return None

        # Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json["metadata"]["sample_name"], f)
        self.add_data_source(f, s_name)

        metrics_dict = parsed_json["metrics"]

        # Add all in the main data_table
        self.mtnuc_data[s_name] = metrics_dict

        # Add version info
        version = parsed_json["metadata"]["version"]
        self.add_software_version(version, s_name)

    def mtnucratio_general_stats_table(self):
        """Take the parsed stats from the mtnucratio report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "mt_cov_avg": {
                "title": "MT genome coverage",
                "description": "Average coverage (X) on mitochondrial genome.",
                "min": 0,
                "scale": "OrRd",
                "suffix": " X",
                "hidden": True,
            },
            "nuc_cov_avg": {
                "title": "Genome coverage",
                "description": "Average coverage (X) on nuclear genome.",
                "min": 0,
                "scale": "GnBu",
                "suffix": " X",
                "hidden": True,
            },
            "mt_nuc_ratio": {
                "title": "MT to Nuclear Ratio",
                "description": "Mitochondrial to nuclear reads ratio (MTNUC)",
                "min": 0,
                "max": 100,
                "scale": "RdYlGn-rev",
            },
            "nucreads": {
                "title": f"{config.read_count_prefix} Genome reads",
                "description": f"Reads on the nuclear genome ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "BuPu",
                "hidden": True,
            },
            "mtreads": {
                "title": f"{config.read_count_prefix} MT genome reads",
                "description": f"Reads on the mitochondrial genome ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "scale": "OrRd",
                "hidden": True,
            },
        }

        self.general_stats_addcols(self.mtnuc_data, headers)
