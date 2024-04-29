""" MultiQC module to parse output from Anglerfish """

import json
import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Anglerfish module class
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Anglerfish",
            anchor="Anglerfish",
            href="https://github.com/remiolsen/anglerfish",
            info="A tool to assess Illumina libraries sequenced on Oxford Nanopore for the purpose of quality control.",
            # doi="", No DOI available
        )

        # Find and load any anglerfish reports
        self.anglerfish_data = dict()

        for f in self.find_log_files("anglerfish", filehandles=True):
            self.parse_anglerfish_json(f)

        # Filter to strip out ignored sample names
        self.anglerfish_data = self.ignore_samples(self.anglerfish_data)

        # Stop execution of the data if no anglerfish data is found.
        if len(self.anglerfish_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.anglerfish_data)} reports")

        # Write parsed report data to a file
        # Parse whole JSON to save all its content
        self.write_data_file(self.anglerfish_data, "multiqc_anglerfish")

        # General Stats Table
        self.anglerfish_general_stats_table()

        # Adds section for Sample Stats Read length table/violin plot
        self.anglerfish_sample_stats()
        # Adds section for Undetermined indexes plot
        self.anglerfish_undetermined_index_chart()

    def parse_anglerfish_json(self, f):
        """Parse the JSON output from Anglerfish and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except Exception:
            file = f["fn"]
            log.warning(f"Could not parse Anglerfish JSON: '{file}'")
            return

        # Fetch a sample name from the command
        s_name = f["s_name"]
        self.add_data_source(f, s_name)
        self.anglerfish_data[s_name] = {}

        # Add version info
        self.add_software_version(parsed_json["anglerfish_version"], s_name)

        # Parse Sample Stats
        ## Index for each sample and their reads in order to iterate without knowing sample names
        index = 0
        total_reads = 0
        try:
            for k in parsed_json["sample_stats"]:
                total_reads += float(k.get("#reads", 0))
                for key in k:
                    try:
                        self.anglerfish_data[s_name][f"{key}_{index}"] = float(k[key])
                    except ValueError:
                        self.anglerfish_data[s_name][f"{key}_{index}"] = k[key]
                index += 1
            self.anglerfish_data[s_name]["sample_stats_amount"] = index
            self.anglerfish_data[s_name]["total_read"] = total_reads
        except KeyError:
            # No sample stat in file or sample stat missing info
            self.anglerfish_data[s_name]["sample_stats_amount"] = -1

        # Parse Undetermined Indexes
        ## Index for each undetermined (count, index) pair
        index = 0
        total_count = 0
        try:
            for k in parsed_json["undetermined"]:
                if len(k) > 0:
                    self.anglerfish_data[s_name][f"undetermined_count_{index}"] = float(k["count"])
                    total_count += float(k["count"])
                    self.anglerfish_data[s_name][f"undetermined_index_{index}"] = k["undetermined_index"]
                    index += 1
            self.anglerfish_data[s_name]["undetermined_amount"] = index
        except KeyError:
            # No undetermined in file or undetermined missing info
            self.anglerfish_data[s_name]["undetermined_amount"] = -1
        self.anglerfish_data[s_name]["total_count"] = total_count

    # General stats table
    def anglerfish_general_stats_table(self):
        """Add Anglerfish statistics to the general statistics table"""
        # Prep data for general stat table
        # Multiple sample names per file requires dict where the first key is not file name
        data = {}
        for s_name in self.anglerfish_data:
            total_read = self.anglerfish_data[s_name]["total_read"]
            total_count = self.anglerfish_data[s_name]["total_count"]
            try:
                for k in range(self.anglerfish_data[s_name]["sample_stats_amount"]):
                    key = self.anglerfish_data[s_name][f"sample_name_{k}"]
                    data[key] = {}
                    data["undetermined"] = {}
                    data[f"total_read_{s_name}"] = {}
                    reads = self.anglerfish_data[s_name][f"#reads_{k}"]
                    data[key]["#reads"] = reads
                    # data[f"total_read_{s_name}"]["#reads"] = total_read
                    data[key]["mean_read_len"] = self.anglerfish_data[s_name][f"mean_read_len_{k}"]
                    data[key]["std_read_len"] = self.anglerfish_data[s_name][f"std_read_len_{k}"]
                    data[key]["library"] = float((reads / total_read) * 100)
                    data["undetermined"]["library"] = float((total_count / total_read) * 100)
            except KeyError:
                log.debug(f"No general stats table generated from Anglerfish json: {s_name}")
            except ZeroDivisionError:
                log.debug(f"No library in general stats table generated from Anglerfish json: {s_name}")

        headers = {
            "library": {
                "title": "% Library",
                "description": "Fraction within library.",
                "max": 100,
                "min": 0,
                "scale": "PuBu-rev",
                "suffix": " %",
            },
            "#reads": {
                "title": "# Reads",
                "description": "Total number of reads",
                "min": 0,
                "scale": "PuOr",
                "format": "{:.0f}",
            },
            "mean_read_len": {
                "title": "Read Length",
                "description": "Mean read length",
                "min": 0,
                "scale": "RdYlGn",
                "suffix": " bp",
            },
            "std_read_len": {
                "title": "Length StdDev",
                "description": "Standard deviation of the read lengths",
                "min": 0,
                "scale": "RdPu",
                "suffix": " bp",
            },
        }

        self.general_stats_addcols(data, headers)

    def anglerfish_sample_stats(self):
        """Generate plot for read length from sample stats.
        For < 10 samples: generate a table
        for >= 10 samples: generate a violin plot"""
        data = {}
        total_samples = 0
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["sample_stats_amount"]
            if index > 0:
                total_samples += index
                for i in range(index):
                    sample_name = self.anglerfish_data[s_name][f"sample_name_{i}"]
                    data[f"Sample: {sample_name}"] = {}
                    data[f"Sample: {sample_name}"]["Mean"] = self.anglerfish_data[s_name][f"mean_read_len_{i}"]
                    data[f"Sample: {sample_name}"]["Standard Deviation"] = self.anglerfish_data[s_name][
                        f"std_read_len_{i}"
                    ]
            else:
                # For non-existing sample stat and faulty sample stat
                log.debug(f"Missing Sample Stat Data in Anglerfish json: {s_name}")
        if len(data) == 0:
            return

        config = {
            "id": "Sample_Stat_Read_Length",
            "title": "Anglerfish: Read Lengths Summary",
        }
        # Plot table if less than 10 samples exist, a violin if more
        if total_samples < 10:
            p = table.plot(data, None, config)
        else:
            p = violin.plot(data, None, config)
        self.add_section(
            name="Read Lengths Summary",
            anchor="anglerfish-sample-statistics",
            description="The Mean read length and the Standard Deviation for each sample.",
            plot=p,
        )

    def anglerfish_undetermined_index_chart(self):
        """Generate Undetermined indexes Bar Plot"""
        data = {}
        for s_name in self.anglerfish_data:
            index = self.anglerfish_data[s_name]["undetermined_amount"]
            # Index smaller than 0 caused by KeyError from no undetermined data
            if index > 0:
                for i in range(index):
                    undetermined_index = self.anglerfish_data[s_name][f"undetermined_index_{i}"]
                    data[undetermined_index] = {}
                    data[undetermined_index][undetermined_index] = self.anglerfish_data[s_name][
                        f"undetermined_count_{i}"
                    ]
            else:
                # For non existing undetermined and faulty undetermined
                log.debug(f"Missing Undetermined Data in Anglerfish json: {s_name}")
        # Only add undetermined section if undetermined data exists
        if len(data) == 0:
            return

        config = {
            "id": "anglerfish_undetermined_index_plot",
            "cpswitch": False,
            "title": "Anglerfish: Undetermined Indexes",
            "ylab": "Index Count",
        }
        self.add_section(
            name="Undetermined Indexes",
            anchor="anglerfish-undetermined-indexes",
            plot=bargraph.plot(data, None, config),
        )
