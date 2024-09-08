import json
import logging
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    AfterQC goes through all FastQ files in a folder and outputs three folders: good, bad and QC folders,
    which contains good reads, bad reads and the QC results of each fastq file/pair.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="AfterQC",
            anchor="afterqc",
            href="https://github.com/OpenGene/AfterQC",
            info="Automatic filtering, trimming, error removing, and quality control for FastQ data.",
            doi="10.1186/s12859-017-1469-3",
        )

        # Find and load any Afterqc reports
        self.afterqc_data: Dict = dict()
        for f in self.find_log_files("afterqc", filehandles=True):
            self.parse_afterqc_log(f)

        # Filter to strip out ignored sample names
        self.afterqc_data = self.ignore_samples(self.afterqc_data)

        if len(self.afterqc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.afterqc_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.afterqc_data, "multiqc_afterqc")

        # Basic Stats Table
        self.afterqc_general_stats_table()

        # Alignment bar plot
        self.add_section(
            name="Bad Reads",
            anchor="after_qc",
            description="Filtering statistics of sampled reads.",
            plot=self.after_qc_bad_reads_chart(),
        )

    def parse_afterqc_log(self, f):
        """Parse the JSON output from AfterQC and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except Exception:
            log.warning(f"Could not parse AfterQC JSON: '{f['fn']}'")
            return None

        # AfterQC changed the name of their summary key at some point
        if "summary" in parsed_json:
            summaryk = "summary"
        elif "afterqc_main_summary" in parsed_json:
            summaryk = "afterqc_main_summary"
        else:
            log.warning(f"AfterQC JSON did not have a 'summary' or 'afterqc_main_summary' key, skipping: '{f['fn']}'")
            return None

        s_name = f["s_name"]
        self.add_data_source(f, s_name)
        self.afterqc_data[s_name] = {}
        for k in parsed_json[summaryk]:
            try:
                self.afterqc_data[s_name][k] = float(parsed_json[summaryk][k])
            except ValueError:
                self.afterqc_data[s_name][k] = parsed_json[summaryk][k]
        try:
            self.afterqc_data[s_name]["pct_good_bases"] = (
                self.afterqc_data[s_name]["good_bases"] / self.afterqc_data[s_name]["total_bases"]
            ) * 100.0
        except KeyError:
            pass

    def afterqc_general_stats_table(self):
        """Take the parsed stats from the Afterqc report and add it to the
        General Statistics table at the top of the report"""

        headers = {
            "pct_good_bases": {
                "title": "% Good Bases",
                "description": "Percent Good Bases",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "BuGn",
            },
            "good_reads": {
                "title": f"{config.read_count_prefix} Good Reads",
                "description": f"Good Reads ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "scale": "GnBu",
                "shared_key": "read_count",
            },
            "total_reads": {
                "title": f"{config.read_count_prefix} Total Reads",
                "description": f"Total Reads ({config.read_count_desc})",
                "min": 0,
                "modify": lambda x: x * config.read_count_multiplier,
                "scale": "Blues",
                "shared_key": "read_count",
            },
            "readlen": {
                "title": "Read Length",
                "description": "Read Length",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.0f}",
                "scale": "YlGn",
            },
        }
        self.general_stats_addcols(self.afterqc_data, headers)

    def after_qc_bad_reads_chart(self):
        """Function to generate the AfterQC bad reads bar plot"""
        # Specify the order of the different possible categories
        keys = {
            "good_reads": {"name": "Good Reads"},
            "bad_reads_with_bad_barcode": {"name": "Bad Barcode"},
            "bad_reads_with_bad_overlap": {"name": "Bad Overlap"},
            "bad_reads_with_bad_read_length": {"name": "Bad Read Length"},
            "bad_reads_with_low_quality": {"name": "Low Quality"},
            "bad_reads_with_polyX": {"name": "PolyX"},
            "bad_reads_with_reads_in_bubble": {"name": "Reads In Bubble"},
            "bad_reads_with_too_many_N": {"name": "Too many N"},
        }

        # Config for the plot
        pconfig = {
            "id": "afterqc_bad_reads_plot",
            "title": "AfterQC: Filtered Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_empty": False,
        }
        return bargraph.plot(self.afterqc_data, keys, pconfig)
