#!/usr/bin/env python

""" MultiQC module to parse output from UMI-tools """


import logging
import os
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    umitools module class, parses dedup logs
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="UMI-tools",
            anchor="umitools",
            href="https://github.com/CGATOxford/UMI-tools",
            info="contains tools for dealing with Unique Molecular Identifiers (UMIs)/(RMTs) and scRNA-Seq barcodes.",
            doi="10.1101/gr.209601.116",
        )

        # Find and load any umitools log files
        self.umitools_data = dict()
        for f in self.find_log_files("umitools/dedup"):
            # Parse the log file for sample name and statistics
            input_fname, data = self.parse_logs(f['f'])
            if len(data) > 1:
                # Clean the sample name
                f["s_name"] = self.clean_s_name(input_fname, f)
                # Log a warning if the log file matches an existing sample name
                if f["s_name"] in self.umitools_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                # Store the data in the overall data dictionary
                self.umitools_data[f["s_name"]] = data
                # Add the log file information to the multiqc_sources.txt
                self.add_data_source(f)

        # Check the log files against the user supplied list of samples to ignore
        self.umitools_data = self.ignore_samples(self.umitools_data)

        # If no files are found, raise an exception
        if len(self.umitools_data) == 0:
            raise UserWarning

        # Log the number of reports found
        log.info("Found {} reports".format(len(self.umitools_data)))

        # Write parsed report data to a file
        self.write_data_file(self.umitools_data, "multiqc_umitools")

        # write data to the general statistics table
        self.umitools_general_stats_table()

        # add a section with a deduplicated reads plot to the report
        self.add_section(
            name="Deduplicated Reads", anchor="umitools-dedup-plot", plot=self.umitools_deduplication_plot()
        )

        # add a section with a table of UMI stats to the report
        self.add_section(name="UMI Stats", anchor="umitools-umi-stats", plot=self.umitools_umi_stats_table())

    def parse_logs(self, f):
        # Initialise a dictionary to hold the data from this log file
        logdata = {}
        # Initialise an empty string for the sample name
        parsed_fname = None

        # Make a dictionary to hold the data lookup table and data types
        lookup_dict = {
            "total_umis": ("INFO total_umis ", int),
            "distinct_umis": ("INFO #umis ", int),
            "input_reads": ("INFO Reads: Input Reads: ", int),
            "output_reads": ("INFO Number of reads out: ", int),
            "positions_deduplicated": ("INFO Total number of positions deduplicated: ", int),
            "mean_umi_per_pos": ("INFO Mean number of unique UMIs per position: ", float),
            "max_umi_per_pos": ("INFO Max. number of unique UMIs per position: ", int),
        }

        # iterate through each line of the log file
        for line in f.splitlines():

            # search for the sample name
            if line.startswith("# stdin"):
                # parse the line and write to the sample name
                parsed_fname = os.path.basename(line.partition("name=")[2].partition(" mode=")[0].strip("'"))

            if parsed_fname is not None:
                # iterate through each item in the lookup table
                for key, value in lookup_dict.items():
                    # search for the lookup value
                    if value[0] in line:
                        # parse the line and write to the data dictionary
                        logdata[key] = value[1](line.partition(value[0])[2])

        # calculate a few simple supplementary stats
        try:
            logdata["percent_passing_dedup"] = round(((logdata["output_reads"] / logdata["input_reads"]) * 100.0), 2)
        except (KeyError, ZeroDivisionError):
            pass
        try:
            logdata["removed_reads"] = logdata["input_reads"] - logdata["output_reads"]
        except (KeyError):
            pass

        return parsed_fname, logdata

    def umitools_general_stats_table(self):
        """Take the parsed stats from the umitools report and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["output_reads"] = {
            "title": "Unique Reads",
            "description": "Reads remaining after deduplication",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "RdYlGn",
        }
        headers["percent_passing_dedup"] = {
            "title": "% Pass Dedup",
            "description": "% processed reads that passed deduplication",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
        }
        self.general_stats_addcols(self.umitools_data, headers)

    def umitools_deduplication_plot(self):
        """Make the HighCharts HTML to plot the deduplication rates"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["output_reads"] = {"color": "#7fc97f", "name": "Reads remaining"}
        keys["removed_reads"] = {"color": "#fdc086", "name": "Reads removed"}

        # Config for the plot
        config = {
            "id": "umitools_deduplication",
            "title": "UMI-tools: Deduplication Counts",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(self.umitools_data, keys, config)

    def umitools_umi_stats_table(self):
        """Take the parsed stats from the umitools report and generate a table of umi stats"""

        headers = OrderedDict()
        headers["positions_deduplicated"] = {
            "title": "Pos Dedup",
            "description": "genomic positions deduplicated",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Greens",
        }
        headers["total_umis"] = {
            "title": "Total UMIs",
            "description": "total umis found in sample",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Blues",
        }
        headers["distinct_umis"] = {
            "title": "Distinct UMIs",
            "description": "distinct umis found in sample",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Purples",
        }
        headers["mean_umi_per_pos"] = {
            "title": "mean #UMI",
            "description": "mean UMIs at each genomic position",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Reds",
        }
        headers["max_umi_per_pos"] = {
            "title": "max #UMI",
            "description": "max UMIs at any genomic position",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Oranges",
        }

        # Config for the table
        config = {
            "id": "umitools_stats_table",
            "table_title": "UMI-tools: UMI stats",
        }

        return table.plot(self.umitools_data, headers, config)
