import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HiCUP",
            anchor="hicup",
            href="http://www.bioinformatics.babraham.ac.uk/projects/hicup/",
            info="Mapping and quality control on Hi-C data.",
            doi="10.12688/f1000research.7334.1",
        )

        # Find and load any HiCUP summary reports
        self.hicup_data = dict()
        for f in self.find_log_files("hicup"):
            self.parse_hicup_logs(f)

        # Filter to strip out ignored sample names
        self.hicup_data = self.ignore_samples(self.hicup_data)

        if len(self.hicup_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hicup_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed data to a file
        self.write_data_file(self.hicup_data, "multiqc_hicup")

        # Basic Stats Table
        self.hicup_stats_table()

        # Report sections
        self.add_section(name="Read Truncation", anchor="hicup-truncating", plot=self.hicup_truncating_chart())
        self.add_section(name="Read Mapping", anchor="hicup-mapping", plot=self.hicup_alignment_chart())

        self.add_section(name="Read Pair Filtering", anchor="hicup-filtering", plot=self.hicup_filtering_chart())

        # TODO: Is there a log file with this data for a line plot?
        # self.add_section (
        #     name = 'Di-Tag Length Distribution',
        #     anchor = 'hicup-lengths',
        #     plot = self.hicup_lengths_chart()
        # )

        self.add_section(
            name="De-Duplication &amp; Di-Tag Separation", anchor="hicup-deduplication", plot=self.hicup_dedup_chart()
        )

    def parse_hicup_logs(self, f):
        """Parse a HiCUP summary report"""
        if not f["fn"].endswith(".txt"):
            return None
        header = []
        lines = f["f"].splitlines()
        for line in lines:
            s = line.split("\t")
            if len(header) == 0:
                if s[0] != "File":
                    return None
                header = s[1:]
            else:
                s_name = self.clean_s_name(s[0], f)
                if s_name.startswith("HiCUP_output/"):
                    s_name = s_name[13:]
                parsed_data = {}
                for idx, num in enumerate(s[1:]):
                    try:
                        parsed_data[header[idx]] = float(num)
                    except ValueError:
                        parsed_data[header[idx]] = num
                parsed_data["Duplicate_Read_Pairs"] = (
                    parsed_data["Valid_Pairs"] - parsed_data["Deduplication_Read_Pairs_Uniques"]
                )
                if s_name in self.hicup_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name)
                self.hicup_data[s_name] = parsed_data

    def hicup_stats_table(self):
        """Add core HiCUP stats to the general stats table"""
        headers = {
            "Percentage_Ditags_Passed_Through_HiCUP": {
                "title": "% Passed",
                "description": "Percentage Di-Tags Passed Through HiCUP",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "Deduplication_Read_Pairs_Uniques": {
                "title": f"{config.read_count_prefix} Unique",
                "description": f"Unique Di-Tags ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "Percentage_Uniques": {
                "title": "% Duplicates",
                "description": "Percent Duplicate Di-Tags",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn-rev",
                "modify": lambda x: 100 - x,
            },
            "Valid_Pairs": {
                "title": f"{config.read_count_prefix} Valid",
                "description": f"Valid Pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "Percentage_Valid": {
                "title": "% Valid",
                "description": "Percent Valid Pairs",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "Paired_Read_1": {
                "title": f"{config.read_count_prefix} Pairs Aligned",
                "description": f"Paired Alignments ({config.read_count_desc})",
                "min": 0,
                "scale": "PuRd",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
            },
            "Percentage_Mapped": {
                "title": "% Aligned",
                "description": "Percentage of Paired Alignments",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
        }
        self.general_stats_addcols(self.hicup_data, headers)

    def hicup_truncating_chart(self):
        """Generate the HiCUP Truncated reads plot"""

        # Specify the order of the different possible categories
        keys = {
            "Not_Truncated_Reads": {"color": "#2f7ed8", "name": "Not Truncated"},
            "Truncated_Read": {"color": "#0d233a", "name": "Truncated"},
        }

        # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
        data = {}
        for s_name in self.hicup_data:
            data[f"{s_name} Read 1"] = {}
            data[f"{s_name} Read 2"] = {}
            data[f"{s_name} Read 1"]["Not_Truncated_Reads"] = self.hicup_data[s_name]["Not_Truncated_Reads_1"]
            data[f"{s_name} Read 2"]["Not_Truncated_Reads"] = self.hicup_data[s_name]["Not_Truncated_Reads_2"]
            data[f"{s_name} Read 1"]["Truncated_Read"] = self.hicup_data[s_name]["Truncated_Read_1"]
            data[f"{s_name} Read 2"]["Truncated_Read"] = self.hicup_data[s_name]["Truncated_Read_2"]

        # Config for the plot
        config = {
            "id": "hicup_truncated_reads_plot",
            "title": "HiCUP: Truncated Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(data, keys, config)

    def hicup_alignment_chart(self):
        """Generate the HiCUP Aligned reads plot"""

        # Specify the order of the different possible categories
        keys = {
            "Unique_Alignments_Read": {"color": "#2f7ed8", "name": "Unique Alignments"},
            "Multiple_Alignments_Read": {"color": "#492970", "name": "Multiple Alignments"},
            "Failed_To_Align_Read": {"color": "#0d233a", "name": "Failed To Align"},
            "Too_Short_To_Map_Read": {"color": "#f28f43", "name": "Too short to map"},
        }

        # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
        data = {}
        for s_name in self.hicup_data:
            data[f"{s_name} Read 1"] = {}
            data[f"{s_name} Read 2"] = {}
            data[f"{s_name} Read 1"]["Unique_Alignments_Read"] = self.hicup_data[s_name]["Unique_Alignments_Read_1"]
            data[f"{s_name} Read 2"]["Unique_Alignments_Read"] = self.hicup_data[s_name]["Unique_Alignments_Read_2"]
            data[f"{s_name} Read 1"]["Multiple_Alignments_Read"] = self.hicup_data[s_name]["Multiple_Alignments_Read_1"]
            data[f"{s_name} Read 2"]["Multiple_Alignments_Read"] = self.hicup_data[s_name]["Multiple_Alignments_Read_2"]
            data[f"{s_name} Read 1"]["Failed_To_Align_Read"] = self.hicup_data[s_name]["Failed_To_Align_Read_1"]
            data[f"{s_name} Read 2"]["Failed_To_Align_Read"] = self.hicup_data[s_name]["Failed_To_Align_Read_2"]
            data[f"{s_name} Read 1"]["Too_Short_To_Map_Read"] = self.hicup_data[s_name]["Too_Short_To_Map_Read_1"]
            data[f"{s_name} Read 2"]["Too_Short_To_Map_Read"] = self.hicup_data[s_name]["Too_Short_To_Map_Read_2"]

        # Config for the plot
        config = {
            "id": "hicup_mapping_stats_plot",
            "title": "HiCUP: Mapping Statistics",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(data, keys, config)

    def hicup_filtering_chart(self):
        """Generate the HiCUP filtering plot"""

        # Specify the order of the different possible categories
        keys = {
            "Valid_Pairs": {"color": "#2f7ed8", "name": "Valid Pairs"},
            "Same_Fragment_Internal": {"color": "#0d233a", "name": "Same Fragment - Internal"},
            "Same_Circularised": {"color": "#910000", "name": "Same Fragment - Circularised"},
            "Same_Dangling_Ends": {"color": "#8bbc21", "name": "Same Fragment - Dangling Ends"},
            "Re_Ligation": {"color": "#1aadce", "name": "Re-ligation"},
            "Contiguous_Sequence": {"color": "#f28f43", "name": "Contiguous Sequence"},
            "Wrong_Size": {"color": "#492970", "name": "Wrong Size"},
        }

        # Config for the plot
        config = {
            "id": "hicup_filtering_plot",
            "title": "HiCUP: Filtering Statistics",
            "ylab": "# Read Pairs",
            "cpswitch_counts_label": "Number of Read Pairs",
            "cpswitch_c_active": False,
        }

        return bargraph.plot(self.hicup_data, keys, config)

    def hicup_dedup_chart(self):
        """Generate the HiCUP Deduplication plot"""

        # Specify the order of the different possible categories
        keys = {
            "Deduplication_Cis_Close_Uniques": {"color": "#2f7ed8", "name": "Unique: cis < 10Kbp"},
            "Deduplication_Cis_Far_Uniques": {"color": "#0d233a", "name": "Unique: cis > 10Kbp"},
            "Deduplication_Trans_Uniques": {"color": "#492970", "name": "Unique: trans"},
            "Duplicate_Read_Pairs": {"color": "#f28f43", "name": "Duplicate read pairs"},
        }

        # Config for the plot
        config = {
            "id": "hicup_dedup_plot",
            "title": "HiCUP: De-Duplication Statistics",
            "ylab": "# Di-Tags",
            "cpswitch_counts_label": "Number of Di-Tags",
            "cpswitch_c_active": False,
        }

        return bargraph.plot(self.hicup_data, keys, config)
