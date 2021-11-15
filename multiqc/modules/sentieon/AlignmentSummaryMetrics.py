#!/usr/bin/env python

""" MultiQC submodule to parse output from Sentieon AlignmentSummaryMetrics
 (based on the Picard module of the same name) """

from collections import OrderedDict
import logging
import os

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Sentieon AlignmentSummaryMetrics reports and parse their data"""

    # Set up vars
    self.sentieon_alignment_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("sentieon/alignment_metrics", filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        for l in f["f"]:
            # New log starting
            if s_name is None and "AlignmentStat" in l:
                keys = None
                # Pull sample name from filename
                s_name = os.path.basename(f["s_name"])
                s_name = self.clean_s_name(s_name, f)
                parsed_data[s_name] = dict()

            if s_name is not None:
                if "AlignmentStat" in l and "#SentieonCommandLine" in l:
                    keys = f["f"].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        # Ignore the FIRST_OF_PAIR / SECOND_OF_PAIR data
                        # to simplify things
                        if vals[0] == "PAIR" or vals[0] == "UNPAIRED":
                            for i, k in enumerate(keys):
                                try:
                                    parsed_data[s_name][k] = float(vals[i])
                                except ValueError:
                                    parsed_data[s_name][k] = vals[i]
                    else:
                        s_name = None
                        keys = None

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # Manipulate sample names if multiple baits found
        for s_name in parsed_data.keys():
            if s_name in self.sentieon_alignment_metrics:
                log.debug(
                    "Duplicate sample name found in {}!\
                          Overwriting: {}".format(
                        f["fn"], s_name
                    )
                )
            self.add_data_source(f, s_name, section="AlignmentSummaryMetrics")
            self.sentieon_alignment_metrics[s_name] = parsed_data[s_name]

    # Filter to strip out ignored sample names
    self.sentieon_alignment_metrics = self.ignore_samples(self.sentieon_alignment_metrics)

    if len(self.sentieon_alignment_metrics) > 0:

        # Write parsed data to a file
        self.write_data_file(self.sentieon_alignment_metrics, "multiqc_sentieon_AlignmentSummaryMetrics")

        # Add to general stats table
        self.general_stats_headers["PCT_PF_READS_ALIGNED"] = {
            "title": "% Aligned",
            "description": "Percent of aligned reads",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.0f}",
            "scale": "RdYlGn",
            "modify": lambda x: self.multiply_hundred(x),
        }
        for s_name in self.sentieon_alignment_metrics:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update(self.sentieon_alignment_metrics[s_name])

        # Make the bar plot of alignment read count
        pdata = dict()
        for s_name in self.sentieon_alignment_metrics.keys():
            pdata[s_name] = dict()
            # Sentieon reports both reads for PE data.
            # Divide it by two as most people will expect # clusters
            if self.sentieon_alignment_metrics[s_name]["CATEGORY"] == "PAIR":
                pdata[s_name]["total_reads"] = self.sentieon_alignment_metrics[s_name]["TOTAL_READS"] / 2
                pdata[s_name]["aligned_reads"] = self.sentieon_alignment_metrics[s_name]["PF_READS_ALIGNED"] / 2
            else:
                pdata[s_name]["total_reads"] = self.sentieon_alignment_metrics[s_name]["TOTAL_READS"]
                pdata[s_name]["aligned_reads"] = self.sentieon_alignment_metrics[s_name]["PF_READS_ALIGNED"]
                pdata[s_name]["unaligned_reads"] = pdata[s_name]["total_reads"] - pdata[s_name]["aligned_reads"]

        keys = OrderedDict()
        keys["aligned_reads"] = {"name": "Aligned Reads"}
        keys["unaligned_reads"] = {"name": "Unaligned Reads"}

        # Config for the plot
        pconfig = {
            "id": "sentieon_aligned_reads",
            "title": "Sentieon: Aligned Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(
            name="Alignment Summary",
            anchor="sentieon-alignmentsummary",
            description="Please note that Sentieon's read counts are divided \
                 by two for paired-end data.",
            plot=bargraph.plot(pdata, keys, pconfig),
        )

    # Return the number of detected samples to the parent module
    return len(self.sentieon_alignment_metrics)
