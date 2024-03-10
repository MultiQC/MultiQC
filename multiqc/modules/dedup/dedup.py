""" MultiQC module to parse output from DeDup """

import json
import logging
from json import JSONDecodeError

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """DeDup module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="DeDup",
            anchor="dedup",
            href="http://www.github.com/apeltzer/DeDup",
            info="is a tool for duplicate removal for merged/collapsed reads in ancient DNA analysis.",
            doi="10.1186/s13059-016-0918-z",
        )

        # Find and load any DeDup reports
        self.dedup_data = dict()

        for f in self.find_log_files("dedup", filehandles=True):
            try:
                self.parseJSON(f)
            except KeyError:
                logging.warning(f"Error loading file {f['fn']}")

        # Filter to strip out ignored sample names
        self.dedup_data = self.ignore_samples(self.dedup_data)

        if len(self.dedup_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.dedup_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.dedup_data, "multiqc_dedup")

        # Basic Stats Table
        self.dedup_general_stats_table()

        # Alignment Rate Plot
        self.add_section(
            description="This plot shows read categories that were either not removed (unique reads) or removed (duplicates).",
            plot=self.dedup_alignment_plot(),
        )

    # Parse our nice little JSON file
    def parseJSON(self, f):
        """Parse the JSON output from DeDup and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
            # Check for Keys existing
            if "metrics" not in parsed_json or "metadata" not in parsed_json:
                log.debug(f"DeDup JSON missing essential keys - skipping sample: '{f['fn']}'")
                return None
        except JSONDecodeError as e:
            log.debug(f"Could not parse DeDup JSON: '{f['fn']}'")
            log.debug(e)
            return None

        # Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json["metadata"]["sample_name"], f)
        self.add_data_source(f, s_name)

        # Add version info
        version = parsed_json["metadata"]["version"]
        self.add_software_version(version, s_name)

        metrics_dict = parsed_json["metrics"]

        for k in metrics_dict:
            metrics_dict[k] = float(metrics_dict[k])

        # Compute (not) removed _mapped_ reads from given values as dedup only affects mapped reads
        # Keep legacy behaviour in case "mapped_reads" cannot be found for <= v0.12.6
        if "mapped_reads" in metrics_dict:
            metrics_dict["mapped_after_dedup"] = (
                metrics_dict["mapped_reads"]
                - metrics_dict["reverse_removed"]
                - metrics_dict["forward_removed"]
                - metrics_dict["merged_removed"]
            )
        else:
            metrics_dict["not_removed"] = (
                metrics_dict["total_reads"]
                - metrics_dict["reverse_removed"]
                - metrics_dict["forward_removed"]
                - metrics_dict["merged_removed"]
            )
        metrics_dict["reads_removed"] = (
            metrics_dict["reverse_removed"] + metrics_dict["forward_removed"] + metrics_dict["merged_removed"]
        )

        # Add all in the main data_table
        self.dedup_data[s_name] = metrics_dict

    def dedup_general_stats_table(self):
        """Take the parsed stats from the DeDup report and add it to the
        basic stats table at the top of the report"""

        ancient_read_count_prefix = getattr(config, "ancient_read_count_prefix", "K")
        ancient_read_count_desc = getattr(config, "ancient_read_count_desc", "thousands")
        ancient_read_count_multiplier = getattr(config, "ancient_read_count_multiplier", 0.001)

        headers = {
            "dup_rate": {
                "title": "Duplication Rate",
                "description": "Percentage of reads categorised as a technical duplicate",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "format": "{:,.0f}",
                "modify": lambda x: x * 100.0,
            },
            "clusterfactor": {
                "title": "ClusterFactor",
                "description": "CF~1 means high library complexity. Large CF means not worth sequencing deeper.",
                "min": 1,
                "max": 100,
                "scale": "OrRd",
                "format": "{:,.2f}",
            },
            "reads_removed": {
                "title": f"{ancient_read_count_prefix} Reads Removed",
                "description": f"Non-unique reads removed after deduplication ({ancient_read_count_desc})",
                "modify": lambda x: x * ancient_read_count_multiplier,
                "shared_key": "read_count",
                "min": 0,
                "hidden": True,
            },
            "mapped_after_dedup": {
                "title": f"{ancient_read_count_prefix} Post-DeDup Mapped Reads",
                "description": f"Unique mapping reads after deduplication ({ancient_read_count_desc})",
                "modify": lambda x: x * ancient_read_count_multiplier,
                "shared_key": "read_count",
                "min": 0,
            },
        }
        self.general_stats_addcols(self.dedup_data, headers)

    def dedup_alignment_plot(self):
        # Specify the order of the different possible categories
        keys = {
            "mapped_after_dedup": {"name": "Unique Retained"},
            "not_removed": {"name": "Not Removed"},
            "reverse_removed": {"name": "Reverse Removed"},
            "forward_removed": {"name": "Forward Removed"},
            "merged_removed": {"name": "Merged Removed"},
        }

        # Config for the plot
        config = {
            "id": "dedup_rates",
            "title": "DeDup: Deduplicated Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.dedup_data, keys, config)
