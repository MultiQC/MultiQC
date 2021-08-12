#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
from collections import OrderedDict
from distutils.version import StrictVersion
import logging
import re

from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Cutadapt module class, parses stderr logs.
    Also understands logs saved by Trim Galore!
    (which contain cutadapt logs)
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Cutadapt",
            anchor="cutadapt",
            href="https://cutadapt.readthedocs.io/",
            info="is a tool to find and remove adapter sequences, primers, poly-A"
            "tails and other types of unwanted sequence from your high-throughput"
            " sequencing reads.",
        )

        # Find and load any Cutadapt reports
        self.cutadapt_data = dict()
        self.cutadapt_length_counts = {"default": dict()}
        self.cutadapt_length_exp = {"default": dict()}
        self.cutadapt_length_obsexp = {"default": dict()}
        self.ends = ["default"]

        for f in self.find_log_files("cutadapt", filehandles=True):
            self.parse_cutadapt_logs(f)

        # Transform trimmed length data by type
        self.transform_trimming_length_data_for_plot()

        # Filter to strip out ignored sample names
        self.cutadapt_data = self.ignore_samples(self.cutadapt_data)

        if len(self.cutadapt_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.cutadapt_data)))

        # Write parsed report data to a file
        self.write_data_file(self.cutadapt_data, "multiqc_cutadapt")

        # Basic Stats Table
        self.cutadapt_general_stats_table()

        # Bar plot with number of reads trimmed
        self.cutadapt_filtered_barplot()

        # Trimming Length Profiles
        self.cutadapt_length_trimmed_plot()

    def parse_cutadapt_logs(self, f):
        """Go through log file looking for cutadapt output"""
        fh = f["f"]
        regexes = {
            "1.7": {
                "bp_processed": "Total basepairs processed:\s*([\d,]+) bp",
                "bp_written": "Total written \(filtered\):\s*([\d,]+) bp",
                "quality_trimmed": "Quality-trimmed:\s*([\d,]+) bp",
                "r_processed": "Total reads processed:\s*([\d,]+)",
                "pairs_processed": "Total read pairs processed:\s*([\d,]+)",
                "r_with_adapters": "Reads with adapters:\s*([\d,]+)",
                "r1_with_adapters": "Read 1 with adapter:\s*([\d,]+)",
                "r2_with_adapters": "Read 2 with adapter:\s*([\d,]+)",
                "r_too_short": "Reads that were too short:\s*([\d,]+)",
                "pairs_too_short": "Pairs that were too short:\s*([\d,]+)",
                "r_too_long": "Reads that were too long:\s*([\d,]+)",
                "pairs_too_long": "Pairs that were too long:\s*([\d,]+)",
                "r_too_many_N": "Reads with too many N:\s*([\d,]+)",
                "pairs_too_many_N": "Pairs with too many N:\s*([\d,]+)",
                "r_written": "Reads written \(passing filters\):\s*([\d,]+)",
                "pairs_written": "Pairs written \(passing filters\):\s*([\d,]+)",
            },
            "1.6": {
                "r_processed": "Processed reads:\s*([\d,]+)",
                "bp_processed": "Processed bases:\s*([\d,]+) bp",
                "r_trimmed": "Trimmed reads:\s*([\d,]+)",
                "quality_trimmed": "Quality-trimmed:\s*([\d,]+) bp",
                "bp_trimmed": "Trimmed bases:\s*([\d,]+) bp",
                "too_short": "Too short reads:\s*([\d,]+)",
                "too_long": "Too long reads:\s*([\d,]+)",
            },
        }
        s_name = None
        end = "default"
        cutadapt_version = None
        parsing_version = "1.7"
        log_section = None
        for l in fh:
            # New log starting
            if "This is cutadapt" in l or "cutadapt version" in l:
                s_name = None
                end = "default"
                cutadapt_version = None
                c_version = re.match(r"This is cutadapt ([\d\.]+)", l)
                if c_version:
                    cutadapt_version = c_version.group(1)
                    try:
                        assert StrictVersion(c_version.group(1)) <= StrictVersion("1.6")
                        parsing_version = "1.6"
                    except:
                        parsing_version = "1.7"
                c_version_old = re.match(r"cutadapt version ([\d\.]+)", l)
                if c_version_old:
                    cutadapt_version = c_version_old.group(1)
                    # The pattern "cutadapt version XX" is only pre-1.6
                    parsing_version = "1.6"
            # Get sample name from end of command line params
            if l.startswith("Command line parameters"):
                for cli in reversed(l.split()):
                    if not cli.startswith("-"):
                        s_name = cli
                        break
                # Manage case where sample name is '-' (reading from stdin)
                if s_name == "-":
                    s_name = f["s_name"]
                else:
                    s_name = self.clean_s_name(s_name, f)
                if s_name in self.cutadapt_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.cutadapt_data[s_name] = dict()
                if cutadapt_version:
                    self.cutadapt_data[s_name]["cutadapt_version"] = cutadapt_version

            if s_name is not None:
                self.add_data_source(f, s_name)

                # Search regexes for overview stats
                for k, r in regexes[parsing_version].items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(",", ""))

                # Starting a new section
                if "===" in l:
                    log_section = l.strip().strip("=").strip()

                # Detect whether 3' or 5'
                end_regex = re.search("Type: regular (\d)'", l)
                if end_regex:
                    end = end_regex.group(1)

                if "Overview of removed sequences" in l:
                    if "' end" in l:
                        res = re.search("(\d)' end", l)
                        end = res.group(1)

                    # Initilise dictionaries for length data if not already done
                    if end not in self.cutadapt_length_counts:
                        self.cutadapt_length_counts[end] = dict()
                        self.cutadapt_length_exp[end] = dict()
                        self.cutadapt_length_obsexp[end] = dict()

                # Histogram showing lengths trimmed
                if "length" in l and "count" in l and "expect" in l:
                    plot_sname = s_name
                    if log_section is not None:
                        plot_sname = "{} - {}".format(s_name, log_section)
                    self.cutadapt_length_counts[end][plot_sname] = dict()
                    self.cutadapt_length_exp[end][plot_sname] = dict()
                    self.cutadapt_length_obsexp[end][plot_sname] = dict()

                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.cutadapt_length_counts[end][plot_sname][a_len] = int(r_seqs.group(2))
                            self.cutadapt_length_exp[end][plot_sname][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.cutadapt_length_obsexp[end][plot_sname][a_len] = float(r_seqs.group(2)) / float(
                                    r_seqs.group(3)
                                )
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.cutadapt_length_obsexp[end][plot_sname][a_len] = float(r_seqs.group(2))
                        else:
                            break
        # Calculate a few extra numbers of our own
        for s_name, d in self.cutadapt_data.items():
            # Percent trimmed
            if "bp_processed" in d and "bp_written" in d:
                self.cutadapt_data[s_name]["percent_trimmed"] = (
                    float(d["bp_processed"] - d["bp_written"]) / d["bp_processed"]
                ) * 100
            elif "bp_processed" in d and "bp_trimmed" in d:
                self.cutadapt_data[s_name]["percent_trimmed"] = (
                    (float(d.get("bp_trimmed", 0)) + float(d.get("quality_trimmed", 0))) / d["bp_processed"]
                ) * 100
            # Add missing filtering categories for pre-1.7 logs
            if StrictVersion(d["cutadapt_version"]) > StrictVersion("1.6"):
                if "r_processed" in d:
                    r_filtered_unexplained = (
                        d["r_processed"]
                        - d.get("r_too_short", 0)
                        - d.get("r_too_long", 0)
                        - d.get("r_too_many_N", 0)
                        - d.get("r_written", 0)
                    )
                    if r_filtered_unexplained > 0:
                        self.cutadapt_data[s_name]["r_filtered_unexplained"] = r_filtered_unexplained
                if "pairs_processed" in d:
                    pairs_filtered_unexplained = (
                        d["pairs_processed"]
                        - d.get("pairs_too_short", 0)
                        - d.get("pairs_too_long", 0)
                        - d.get("pairs_too_many_N", 0)
                        - d.get("pairs_written", 0)
                    )
                    if pairs_filtered_unexplained > 0:
                        self.cutadapt_data[s_name]["pairs_filtered_unexplained"] = pairs_filtered_unexplained

    def transform_trimming_length_data_for_plot(self):
        """Check if we parsed double ended data and transform it accordingly"""
        # Sanity check
        if (
            self.cutadapt_length_counts.keys() != self.cutadapt_length_exp.keys()
            or self.cutadapt_length_exp.keys() != self.cutadapt_length_obsexp.keys()
        ):
            log.error("Something went wrong...")
            log.debug("Keys in trimmed length data differed")
            raise UserWarning

        if len(self.cutadapt_length_counts["default"]) == 0:
            self.cutadapt_length_counts.pop("default")
            self.cutadapt_length_exp.pop("default")
            self.cutadapt_length_obsexp.pop("default")

        self.ends = list(self.cutadapt_length_counts.keys())

    def cutadapt_general_stats_table(self):
        """Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report"""

        headers = {}
        headers["percent_trimmed"] = {
            "title": "% BP Trimmed",
            "description": "% Total Base Pairs trimmed",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlBu-rev",
        }
        self.general_stats_addcols(self.cutadapt_data, headers)

    def cutadapt_filtered_barplot(self):
        """Bar plot showing proportion of reads trimmed"""

        pconfig = {"id": "cutadapt_filtered_reads_plot", "title": "Cutadapt: Filtered Reads", "ylab": "Counts"}

        # We just use all categories. If a report is generated with a mixture
        # of SE and PE data then this means quite a lot of categories.
        # Usually, only a single data type is used though - in that case
        # any categories with 0 across all samples will be ignored.
        cats = OrderedDict()
        cats["pairs_written"] = {"name": "Pairs passing filters"}
        cats["r_written"] = {"name": "Reads passing filters"}
        cats["pairs_too_short"] = {"name": "Pairs that were too short"}
        cats["r_too_short"] = {"name": "Reads that were too short"}
        cats["pairs_too_long"] = {"name": "Pairs that were too long"}
        cats["r_too_long"] = {"name": "Reads that were too long"}
        cats["pairs_too_many_N"] = {"name": "Pairs with too many N"}
        cats["r_too_many_N"] = {"name": "Reads with too many N"}
        cats["pairs_filtered_unexplained"] = {"name": "Filtered pairs (uncategorised)"}
        cats["r_filtered_unexplained"] = {"name": "Filtered reads (uncategorised)"}

        self.add_section(
            name="Filtered Reads",
            anchor="cutadapt_filtered_reads",
            description="This plot shows the number of reads (SE) / pairs (PE) removed by Cutadapt.",
            plot=bargraph.plot(self.cutadapt_data, cats, pconfig),
        )

    def cutadapt_length_trimmed_plot(self):
        """Generate the trimming length plot"""
        for end in [x for x in ["default", "5", "3"] if x in self.ends]:
            pconfig = {
                "id": f"cutadapt_trimmed_sequences_plot_{end}",
                "title": "Cutadapt: Lengths of Trimmed Sequences{}".format(
                    "" if end == "default" else f" ({end}' end)"
                ),
                "ylab": "Counts",
                "xlab": "Length Trimmed (bp)",
                "xDecimals": False,
                "ymin": 0,
                "tt_label": "<b>{point.x} bp trimmed</b>: {point.y:.0f}",
                "data_labels": [
                    {"name": "Counts", "ylab": "Count"},
                    {"name": "Obs/Exp", "ylab": "Observed / Expected"},
                ],
            }

            self.add_section(
                name="Trimmed Sequence Lengths{}".format("" if end == "default" else f" ({end}')"),
                anchor="cutadapt_trimmed_sequences{}".format("" if end == "default" else f"_{end}"),
                description="This plot shows the number of reads with certain lengths of adapter trimmed{}.".format(
                    "" if end == "default" else f" for the {end}' end"
                ),
                helptext="""
                Obs/Exp shows the raw counts divided by the number expected due to sequencing errors.
                A defined peak may be related to adapter length.

                See the [cutadapt documentation](http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report)
                for more information on how these numbers are generated.
                """,
                plot=linegraph.plot([self.cutadapt_length_counts[end], self.cutadapt_length_obsexp[end]], pconfig),
            )
