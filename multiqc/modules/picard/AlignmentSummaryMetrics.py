""" MultiQC submodule to parse output from Picard AlignmentSummaryMetrics """

import logging
import os
import re
from collections import OrderedDict

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard AlignmentSummaryMetrics reports and parse their data"""

    # Set up vars
    self.picard_alignment_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(f"{self.anchor}/alignment_metrics", filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        for l in f["f"]:
            # New log starting
            if "AlignmentSummaryMetrics" in l and "INPUT" in l:
                s_name = None
                keys = None
                # Pull sample name from input
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip("[]"))
                    s_name = self.clean_s_name(s_name, f)

            if s_name is None and "AlignmentStat" in l:
                # Pull sample name from filename
                s_name = os.path.basename(f["s_name"])
                s_name = self.clean_s_name(s_name, f)
                parsed_data[s_name] = dict()
                keys = None

            if ("AlignmentSummaryMetrics" in l and "## METRICS CLASS" in l) or (
                "AlignmentStat" in l and "#SentieonCommandLine" in l
            ):
                keys = f["f"].readline().strip("\n").split("\t")
            elif keys:
                if s_name not in parsed_data:
                    parsed_data[s_name] = dict()
                vals = l.strip("\n").split("\t")
                if len(vals) == len(keys):
                    # Ignore the FIRST_OF_PAIR / SECOND_OF_PAIR data to simplify things
                    if vals[0] == "PAIR" or vals[0] == "UNPAIRED":
                        for i, k in enumerate(keys):
                            try:
                                parsed_data[s_name][k] = float(vals[i])
                            except ValueError:
                                parsed_data[s_name][k] = vals[i]
                else:
                    s_name = None
                    keys = None

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None, s_name)

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # Manipulate sample names if multiple baits found
        for s_name in parsed_data.keys():
            if s_name in self.picard_alignment_metrics:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
            self.add_data_source(f, s_name, section="AlignmentSummaryMetrics")
            self.picard_alignment_metrics[s_name] = parsed_data[s_name]

    # Filter to strip out ignored sample names
    self.picard_alignment_metrics = self.ignore_samples(self.picard_alignment_metrics)

    if len(self.picard_alignment_metrics) > 0:
        # Write parsed data to a file
        self.write_data_file(self.picard_alignment_metrics, f"multiqc_{self.anchor}_AlignmentSummaryMetrics")

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
        for s_name in self.picard_alignment_metrics:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update(self.picard_alignment_metrics[s_name])

        # Make the bar plot of alignment read count + # aligned bases
        pdata = dict()
        for s_name in self.picard_alignment_metrics.keys():
            pdata[s_name] = dict()
            # Picard reports both reads for PE data. Divide it by two as most people will expect # clusters
            if self.picard_alignment_metrics[s_name]["CATEGORY"] == "PAIR":
                pdata[s_name]["total_reads"] = self.picard_alignment_metrics[s_name]["TOTAL_READS"] / 2
                pdata[s_name]["aligned_reads"] = self.picard_alignment_metrics[s_name]["PF_READS_ALIGNED"] / 2
            else:
                pdata[s_name]["total_reads"] = self.picard_alignment_metrics[s_name]["TOTAL_READS"]
                pdata[s_name]["aligned_reads"] = self.picard_alignment_metrics[s_name]["PF_READS_ALIGNED"]
            pdata[s_name]["unaligned_reads"] = pdata[s_name]["total_reads"] - pdata[s_name]["aligned_reads"]

        keys = [OrderedDict(), OrderedDict()]
        keys[0]["aligned_reads"] = {"name": "Aligned Reads"}
        keys[0]["unaligned_reads"] = {"name": "Unaligned Reads"}
        keys[1]["PF_ALIGNED_BASES"] = {"name": "Aligned Bases"}

        # Config for the plot
        pconfig = {
            "id": f"{self.anchor}_alignment_summary",
            "title": f"{self.name}: Alignment Summary",
            "ylab": "# Reads",
            "data_labels": [
                {
                    "name": "Aligned Reads",
                    "ylab": "# Reads",
                    "cpswitch_counts_label": "Number of Reads",
                },
                {
                    "name": "Aligned Bases",
                    "ylab": "# Bases",
                    "cpswitch_counts_label": "Number of Bases",
                },
            ],
        }

        # The different data sets we want to plot
        self.add_section(
            name="Alignment Summary",
            anchor=f"{self.anchor}-alignmentsummary",
            description=f"Please note that {self.name}'s read counts are divided by two for paired-end data. Total bases (including unaligned) is not provided.",
            plot=bargraph.plot([pdata, self.picard_alignment_metrics], keys, pconfig),
        )

        # Make a bar plot of mean read length
        keys = {"MEAN_READ_LENGTH": {"name": "Mean Read Length"}}
        pconfig = {
            "id": f"{self.anchor}_alignment_readlength_plot",
            "title": f"{self.name}: Mean Read Length",
            "ylab": "Base pairs",
            "cpswitch": False,
        }

        # The different data sets we want to plot
        self.add_section(
            name="Mean read length",
            anchor=f"{self.anchor}_alignment_readlength",
            description="The mean read length of the set of reads examined.",
            plot=bargraph.plot(self.picard_alignment_metrics, keys, pconfig),
        )

    # Return the number of detected samples to the parent module
    return len(self.picard_alignment_metrics)
