""" MultiQC module to parse output from pycoQC """

import logging

import yaml

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="pycoQC",
            anchor="pycoqc",
            href="https://github.com/tleonardi/pycoQC",
            info="computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data",
            doi="10.21105/joss.01236",
        )

        self.pycoqc_data = {}
        for f in self.find_log_files("pycoqc"):
            if f["s_name"] in self.pycoqc_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")
            data = self.load_data(f["f"])
            # Function can return None if YAML parsing failed
            if data:
                self.pycoqc_data[f["s_name"]] = data
                self.add_software_version(data["pycoqc"]["version"], f["s_name"])
            self.add_data_source(f)

        self.pycoqc_data = self.ignore_samples(self.pycoqc_data)

        # Stop if we didn't find anything
        if len(self.pycoqc_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.pycoqc_data)} reports")

        # Write data to file
        self.write_data_file(self.pycoqc_data, "pycoqc")

        self.parse_data()

        self.make_general_stats()
        self.make_pycoqc_table()
        self.read_base_count_plot()

        # Only add these sections if we have data (not available for older versions of PycoQC)
        if self.read_length_plot_data[0]:
            self.read_length_plot()

        if self.quality_plot_data[0]:
            self.make_quality_plot()

    def load_data(self, f):
        """Load the PycoQC YAML file"""
        try:
            return yaml.load(f, Loader=yaml.SafeLoader)
        except Exception as e:
            log.warning(f"Could not parse YAML for '{f}': \n  {e}")
            return None

    def parse_data(self):
        """Convert the parsed data into the correct structures for MultiQC functions"""

        self.table_data = {}
        self.reads_data = {}
        self.bases_data = {}
        length_plot_all = {}
        length_plot_pass = {}
        qual_plot_all = {}
        qual_plot_pass = {}

        for sample, sample_data in self.pycoqc_data.items():
            self.table_data[sample] = {
                "all_median_read_length": sample_data["All Reads"]["basecall"]["len_percentiles"][50],
                "all_median_phred_score": sample_data["All Reads"]["basecall"]["qual_score_percentiles"][50],
                "all_n50": sample_data["All Reads"]["basecall"]["N50"],
                "all_run_duration": sample_data["All Reads"]["run"]["run_duration"],
                "all_channels": sample_data["All Reads"]["run"]["active_channels"],
                "all_reads": sample_data["All Reads"]["basecall"]["reads_number"],
                "all_bases": sample_data["All Reads"]["basecall"]["bases_number"],
                "passed_median_read_length": sample_data["Pass Reads"]["basecall"]["len_percentiles"][50],
                "passed_median_phred_score": sample_data["Pass Reads"]["basecall"]["qual_score_percentiles"][50],
                "passed_n50": sample_data["Pass Reads"]["basecall"]["N50"],
                "passed_channels": sample_data["Pass Reads"]["run"]["active_channels"],
                "passed_reads": sample_data["Pass Reads"]["basecall"]["reads_number"],
                "passed_bases": sample_data["Pass Reads"]["basecall"]["bases_number"],
            }

            self.reads_data[sample] = {
                "passed_reads": sample_data["Pass Reads"]["basecall"]["reads_number"],
                "non_passed_reads": sample_data["All Reads"]["basecall"]["reads_number"]
                - sample_data["Pass Reads"]["basecall"]["reads_number"],
            }

            self.bases_data[sample] = {
                "passed_bases": sample_data["Pass Reads"]["basecall"]["bases_number"],
                "non_passed_bases": sample_data["All Reads"]["basecall"]["bases_number"]
                - sample_data["Pass Reads"]["basecall"]["bases_number"],
            }
            try:
                length_x_vals_all = sample_data["All Reads"]["basecall"]["len_hist"]["x"]
                length_y_vals_all = sample_data["All Reads"]["basecall"]["len_hist"]["y"]
                length_plot_all[sample] = dict(zip(length_x_vals_all, length_y_vals_all))
                length_x_vals_pass = sample_data["Pass Reads"]["basecall"]["len_hist"]["x"]
                length_y_vals_pass = sample_data["Pass Reads"]["basecall"]["len_hist"]["y"]
                length_plot_pass[sample] = dict(zip(length_x_vals_pass, length_y_vals_pass))

                qual_x_vals_all = sample_data["All Reads"]["basecall"]["qual_score_hist"]["x"]
                qual_y_vals_all = sample_data["All Reads"]["basecall"]["qual_score_hist"]["y"]
                qual_plot_all[sample] = dict(zip(qual_x_vals_all, qual_y_vals_all))
                qual_x_vals_pass = sample_data["Pass Reads"]["basecall"]["qual_score_hist"]["x"]
                qual_y_vals_pass = sample_data["Pass Reads"]["basecall"]["qual_score_hist"]["y"]
                qual_plot_pass[sample] = dict(zip(qual_x_vals_pass, qual_y_vals_pass))
            except KeyError:
                log.debug(
                    "'{}': Could not find plot data. Please make sure you are using pycoQC v2.5.0.20 or newer.".format(
                        sample
                    )
                )

        self.read_length_plot_data = [length_plot_pass, length_plot_all]
        self.quality_plot_data = [qual_plot_pass, qual_plot_all]

    def make_general_stats(self):
        general_stats_headers = {
            "passed_median_read_length": {
                "title": "Read Length - Pass (bp)",
                "description": "Median Read Length - passing reads (base pairs)",
                "scale": "BuPu",
                "shared_key": "median_read_len",
                "format": "{:,.0f}",
            },
            "all_median_read_length": {
                "title": "Read Length - All (bp)",
                "description": "Median Read Length - all reads (base pairs)",
                "scale": "BuPu",
                "shared_key": "median_read_len",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "passed_reads": {
                "title": f"{config.long_read_count_prefix} Reads - Pass",
                "description": f"Number of reads - passing reads ({config.long_read_count_desc})",
                "scale": "BuGn",
                "modify": lambda x: x * config.long_read_count_multiplier,
                "shared_key": "long_read_count",
            },
            "all_reads": {
                "title": f"{config.long_read_count_prefix} Reads - All",
                "description": f"Number of reads - all reads ({config.long_read_count_desc})",
                "scale": "BuGn",
                "modify": lambda x: x * config.long_read_count_multiplier,
                "shared_key": "long_read_count",
                "hidden": True,
            },
            "passed_bases": {
                "title": f"{config.base_count_prefix} Bases - Pass",
                "description": f"Number of bases - passing reads ({config.base_count_desc} of base pairs)",
                "scale": "OrRd",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
            },
            "all_bases": {
                "title": f"{config.base_count_prefix} Bases - All",
                "description": f"Number of bases - all reads ({config.base_count_desc} of base pairs)",
                "scale": "OrRd",
                "modify": lambda x: x * config.base_count_multiplier,
                "shared_key": "base_count",
                "hidden": True,
            },
        }

        self.general_stats_addcols(self.table_data, general_stats_headers)

    def make_pycoqc_table(self):
        pycoqc_table_headers = {
            "passed_n50": {
                "title": "N50 - Pass (bp)",
                "description": "N50 - passing reads (base pairs)",
                "scale": "BuGn",
                "shared_key": "n50",
                "format": "{:,.0f}",
            },
            "all_n50": {
                "title": "N50 - All (bp)",
                "description": "N50 - all reads (base pairs)",
                "scale": "BuGn",
                "shared_key": "n50",
                "format": "{:,.0f}",
            },
            "passed_median_phred_score": {
                "title": "Median read qual - Pass",
                "description": "Median PHRED quality score - passing reads",
                "scale": "PuRd",
                "shared_key": "phred",
            },
            "all_median_phred_score": {
                "title": "Median read qual - All",
                "description": "Median PHRED quality score - all reads",
                "scale": "PuRd",
                "shared_key": "phred",
            },
            "passed_channels": {
                "title": "Active Channels - Pass",
                "description": "Number of active channels - passing reads",
                "scale": "OrRd",
                "shared_key": "channels",
                "format": "{:,.0f}",
            },
            "all_channels": {
                "title": "Active Channels - All",
                "description": "Number of active channels - all reads",
                "scale": "OrRd",
                "shared_key": "channels",
                "format": "{:,.0f}",
            },
            "all_run_duration": {
                "title": "Run duration (h)",
                "description": "Run duration (hours)",
                "scale": "PuBuGn",
            },
        }

        pycoqc_table_config = {
            "id": "pycoqc_stats_table",
            "namespace": "pycoQC",
            "title": "pycoQC: Statistics",
        }

        self.add_section(
            name="Statistics",
            anchor="pycoqc_stats",
            plot=table.plot(self.table_data, pycoqc_table_headers, pycoqc_table_config),
        )

    def read_base_count_plot(self):
        """Make a section with a bar plot showing the counts of reads and bases
        before and after filtering"""

        # Plot config
        read_bar_config = {
            "id": "pycoqc_count_plot",
            "title": "pycoQC: Read Counts",
            "ylab": "Sample",
            "data_labels": [
                {"name": "Reads", "ylab": "Number of Reads"},
                {"name": "Bases", "ylab": "Number of Base Pairs"},
            ],
        }

        # Two sets of categories for the two datasets
        read_cats = {
            "passed_reads": {"name": "Passed Reads"},
            "non_passed_reads": {"name": "Failed Reads"},
        }
        base_cats = {
            "passed_bases": {"name": "Passed Bases"},
            "non_passed_bases": {"name": "Failed Bases"},
        }

        # Make the plot
        bargraph_plot = bargraph.plot([self.reads_data, self.bases_data], [read_cats, base_cats], read_bar_config)

        self.add_section(
            name="Read / Base counts",
            anchor="pycoqc_counts",
            description="Number of sequenced reads / bases passing and failing QC thresholds.",
            plot=bargraph_plot,
        )

    def read_length_plot(self):
        """Read length plot, showing passed reads and all reads"""

        read_length_config = {
            "id": "pycoqc_read_len_plot",
            "title": "pycoQC: Read Length",
            "ylab": "Read Density",
            "xlab": "Basecalled Length (bp)",
            "xLog": True,
            "data_labels": [
                {"name": "Passing Reads", "ylab": "Read Density", "xlab": "Basecalled Length (bp)"},
                {"name": "All Reads", "ylab": "Read Density", "xlab": "Basecalled Length (bp)"},
            ],
        }
        self.add_section(
            name="Read length",
            anchor="pycoqc_read_len",
            description="Distribution of read length for all / passed reads.",
            plot=linegraph.plot(self.read_length_plot_data, read_length_config),
        )

    def make_quality_plot(self):
        """Quality plot, showing distrubtion of PHRED scores."""

        qual_config = {
            "id": "pycoqc_read_qual_plot",
            "title": "pycoQC: Read Quality Scores",
            "ylab": "Read Density",
            "xlab": "Read Quality PHRED Score",
            "xmin": 0,
            "data_labels": [
                {"name": "Passing Reads", "ylab": "Read Density", "xlab": "Read Quality PHRED Score"},
                {"name": "All Reads", "ylab": "Read Density", "xlab": "Read Quality PHRED Score"},
            ],
        }
        self.add_section(
            name="Quality scores",
            anchor="pycoqc_read_qual",
            description="Distribution of quality scores for all / passed reads.",
            plot=linegraph.plot(self.quality_plot_data, qual_config),
        )
