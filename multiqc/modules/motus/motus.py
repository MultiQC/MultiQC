""" Module to parse output from mOTUs """

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import config

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="motus",
            anchor="motus",
            href="https://motu-tool.org/",
            info="is a tool performing microbial profiling through marker gene (MG)-based operational taxonomic units (mOTUs).",
            doi="10.1038/s41467-019-08844-4",
        )

        # Define the main motus multiqc data object
        self.motus_data = dict()

        for f in self.find_log_files("motus", filehandles=True):
            self.parse_logs(f)

        self.motus_data = self.ignore_samples(self.motus_data)

        if len(self.motus_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.motus_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write data to file
        self.write_data_file(self.motus_data, "motus")

        self.motus_general_stats()
        self.motus_filtering_bargraph_plot()
        self.motus_mapping_bargraph_plot()
        self.motus_motus_bargraph_plot()

    def parse_logs(self, f):
        """Parses a motus stdout saved in a file"""
        # Assume take from the name of the file as not always reported in the log
        s_name = f["s_name"]

        if self.motus_data.get(s_name) is not None:
            log.warning(f"Duplicate sample name found based on filename! Overwriting: {s_name}")

        self.motus_data[s_name] = {}
        self.add_data_source(f, s_name=s_name)
        for line in f["f"]:
            # Search strings - keeping colon as to not pick up non-stats lines
            for cat in [
                "Total number of reads:",
                "Total number of inserts:",
                "Unique mappers:",
                "Multiple mappers:",
                "Ignored multiple mapper without unique hit:",
                "Number of ref-mOTUs:",
                "Number of meta-mOTUs:",
                "Number of ext-mOTUs:",
            ]:
                if cat in line:
                    self.motus_data[s_name][cat.replace(":", "")] = int(line.strip().split(":")[1].lstrip())
                elif "Number of reads after filtering:" in line:
                    self.motus_data[s_name]["Number of reads after filtering"] = int(
                        line.strip().split(":")[1].lstrip().split(" ")[0]
                    )
                    self.motus_data[s_name]["Percent reads after filtering"] = float(
                        line.strip().split(":")[1].lstrip().split(" ")[1].replace("(", "")
                    )
        self.motus_data[s_name]["Discarded reads"] = (
            self.motus_data[s_name]["Total number of reads"]
            - self.motus_data[s_name]["Number of reads after filtering"]
        )

    def motus_general_stats(self):
        """mOTUs read counts for general stats"""
        headers = {
            "Total number of reads": {
                "title": f"Total Input Reads ({config.read_count_prefix})",
                "description": f"Total number of input reads to mOTUs ({config.read_count_prefix})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "Number of reads after filtering": {
                "title": f"Total Mapped Reads ({config.read_count_prefix})",
                "description": f"Total number of reads after mapping({config.read_count_prefix})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "Total number of inserts": {
                "title": f"Total Mapped Inserts ({config.read_count_prefix})",
                "description": f"Total number of inserts mapped to a MGC ({config.read_count_prefix})",
                "scale": "Purples",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "Unique mappers": {
                "title": f"Unique Mapped Inserts ({config.read_count_prefix})",
                "description": f"Total number of inserts mapped to a single MGC ({config.read_count_prefix})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "Multiple mappers": {
                "title": f" Multi-mapped Inserts ({config.read_count_prefix})",
                "description": f"Total number of inserts mapped to multiple MGCs ({config.read_count_prefix})",
                "scale": "Purples",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "Ignored multiple mapper without unique hit": {
                "title": f"Ignored Multi-mapped Inserts ({config.read_count_prefix})",
                "description": f"Total number of ignored multi-MGC mapped reads ({config.read_count_prefix})",
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
            "Number of ref-mOTUs": {
                "title": "Total ref-mOTUs",
                "description": "Total known species mOTUs found",
                "scale": "Purples",
                "shared_key": "mOTUs",
                "format": "{:,.0f}",
            },
            "Number of meta-mOTUs": {
                "title": "Total meta-mOTUs",
                "description": "Total number of unknown metagenome-derived mOTUs found",
                "scale": "Greens",
                "shared_key": "mOTUs",
                "format": "{:,.0f}",
            },
            "Number of ext-mOTUs": {
                "title": "Total ext-mOTUs",
                "description": "Total number of unknown MAG-derived mOTUs found",
                "scale": "Purples",
                "shared_key": "mOTUs",
                "format": "{:,.0f}",
            },
        }

        self.general_stats_addcols(self.motus_data, headers)

    def motus_filtering_bargraph_plot(self):
        """mOTUs read counts for general stats"""
        common = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": f"{config.read_count_prefix} reads",
            "tt_decimals": 0,
            "shared_key": "read_count",
        }
        cats = {
            "Number of reads after filtering": dict(common, **{"name": "Reads after mapping"}),
            "Discarded reads": dict(common, **{"name": "Unmapped reads"}),
        }

        self.add_section(
            name="mOTUs: Read filtering information",
            anchor="motus-filtering",
            description="Read filtering statistics (i.e. mapping of reads to the mOTUs marker database).",
            plot=bargraph.plot(
                self.motus_data,
                cats,
                {
                    "id": "motus-filtering-reads",
                    "title": "mOTUs: Read filtering information",
                    "ylab": "Reads",
                },
            ),
        )

    def motus_mapping_bargraph_plot(self):
        """mOTUs bar chart of insert types"""
        common = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": f"{config.read_count_prefix} reads",
            "tt_decimals": 0,
            "shared_key": "read_count",
        }
        cats = {
            "Unique mappers": dict(common, **{"name": "Unique mapped inserts", "color": "#3aba5e"}),
            "Multiple mappers": dict(common, **{"name": "Multiple mapped inserts", "color": "#ebbe59"}),
            "Ignored multiple mapper without unique hit": dict(
                common, **{"name": "Ignored multi-mapped inserts", "color": "#cf5565"}
            ),
        }

        self.add_section(
            name="mOTUs: Insert mapping information",
            anchor="motus-mapping",
            description="How inserts was classified after alignment to MGCs.",
            plot=bargraph.plot(
                self.motus_data,
                cats,
                {
                    "id": "motus-mapping-inserts",
                    "title": "mOTUs: Insert mapping information",
                    "ylab": "Inserts",
                },
            ),
        )

    def motus_motus_bargraph_plot(self):
        """mOTUs bar chart of mOTU types"""
        common = {
            "min": 0,
            "tt_decimals": 0,
        }
        cats = {
            "Number of ref-mOTUs": dict(common, **{"name": "Known mOTUs"}),
            "Number of meta-mOTUs": dict(common, **{"name": "(Unknown) Metagenome mOTUs"}),
            "Number of ext-mOTUs": dict(common, **{"name": "(Unknown) MAG mOTUs"}),
        }

        self.add_section(
            name="mOTUs: mOTU identification information",
            anchor="motus-identification",
            description="Distribution of the types of mOTUs found.",
            plot=bargraph.plot(
                self.motus_data,
                cats,
                {
                    "id": "motus-identification-types",
                    "title": "mOTUs: mOTU identification information",
                    "ylab": "mOTUs",
                },
            ),
        )
