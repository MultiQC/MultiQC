""" MultiQC submodule to parse output from Picard AlignmentSummaryMetrics """

import logging
from collections import OrderedDict

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """
    Find Picard AlignmentSummaryMetrics reports and parse their data.
    """

    # Set up vars
    self.picard_alignment_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(f"{self.anchor}/alignment_metrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        parsed_data = {s_name: dict()}
        # A file can be concatenated from multiple samples, so we need to keep track of
        # the current sample name and header.
        keys = None

        for l in f["f"]:
            maybe_s_name = self.extract_sample_name(l, f, "AlignmentSummaryMetrics")
            if maybe_s_name:
                # Starts information for a new sample
                s_name = maybe_s_name
                keys = None
                if s_name in self.picard_alignment_metrics or s_name in parsed_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                parsed_data[s_name] = dict()

            if s_name is None:
                continue

            if self.is_line_right_before_table(l):
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

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # When there is only one sample, using the file name to extract the sample name.
        # if len(parsed_data) == 1:
        #     parsed_data = {f["s_name"]: list(parsed_data.values())[0]}

        self.picard_alignment_metrics.update(parsed_data)

        for s_name in parsed_data:
            self.add_data_source(f, s_name, section="AlignmentSummaryMetrics")

    # Filter to strip out ignored sample names
    self.picard_alignment_metrics = self.ignore_samples(self.picard_alignment_metrics)

    if len(self.picard_alignment_metrics) == 0:
        return 0

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
