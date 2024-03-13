""" MultiQC module to parse output from Cutadapt """


import logging
import os
import re
import shlex
import json

from packaging import version

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

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
            info="""is a tool to find and remove adapter sequences, primers, poly-A
                    tails and other types of unwanted sequence from your high-throughput
                    sequencing reads.""",
            doi="10.14806/ej.17.1.200",
        )

        # Find and load any Cutadapt reports
        self.cutadapt_data = dict()
        self.cutadapt_length_counts = {"default": dict()}
        self.cutadapt_length_exp = {"default": dict()}
        self.cutadapt_length_obsexp = {"default": dict()}

        for f in self.find_log_files("cutadapt"):
            self.parse_file(f)

        # Transform trimmed length data by type
        self.transform_trimming_length_data_for_plot()

        # Filter to strip out ignored sample names
        self.cutadapt_data = self.ignore_samples(self.cutadapt_data)

        if len(self.cutadapt_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.cutadapt_data)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.cutadapt_data, "multiqc_cutadapt")

        # Basic Stats Table
        self.cutadapt_general_stats_table()

        # Bar plot with number of reads trimmed
        self.cutadapt_filtered_barplot()

        # Trimming Length Profiles
        self.cutadapt_length_trimmed_plot()

    def parse_file(self, f):
        if f["fn"].endswith(".json"):
            self.parse_json(f)
        else:
            self.parse_log(f)

        # Calculate a few extra numbers of our own
        for s_name, d in self.cutadapt_data.items():
            # Percent trimmed
            if d.get("bp_processed"):
                if d.get("bp_written") is not None:
                    self.cutadapt_data[s_name]["percent_trimmed"] = (
                        float(d["bp_processed"] - d["bp_written"]) / d["bp_processed"]
                    ) * 100
                elif d.get("bp_trimmed") is not None:
                    self.cutadapt_data[s_name]["percent_trimmed"] = (
                        (float(d.get("bp_trimmed", 0)) + float(d.get("quality_trimmed", 0))) / d["bp_processed"]
                    ) * 100
            # Add missing filtering categories for pre-1.7 logs
            if version.parse(d["cutadapt_version"]) > version.parse("1.6"):
                if d.get("r_processed") is not None:
                    r_filtered_unexplained = (
                        d["r_processed"]
                        - (d.get("r_too_short") or 0)
                        - (d.get("r_too_long") or 0)
                        - (d.get("r_too_many_N") or 0)
                        - (d.get("r_written") or 0)
                    )
                    if r_filtered_unexplained > 0:
                        self.cutadapt_data[s_name]["r_filtered_unexplained"] = r_filtered_unexplained
                if d.get("pairs_processed") is not None:
                    pairs_filtered_unexplained = (
                        d["pairs_processed"]
                        - (d.get("pairs_too_short") or 0)
                        - (d.get("pairs_too_long") or 0)
                        - (d.get("pairs_too_many_N") or 0)
                        - (d.get("pairs_written") or 0)
                    )
                    if pairs_filtered_unexplained > 0:
                        self.cutadapt_data[s_name]["pairs_filtered_unexplained"] = pairs_filtered_unexplained

    def parse_json(self, f):
        path = os.path.join(f["root"], f["fn"])
        with open(path, "r") as fh:
            data = json.load(fh)

        s_name = self.clean_s_name([v for k, v in data["input"].items() if k.startswith("path")], f)
        if s_name in self.cutadapt_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

        self.add_software_version(data["cutadapt_version"], s_name)
        self.add_data_source(f, s_name)

        d = dict()
        d["cutadapt_version"] = data["cutadapt_version"]
        d["bp_processed"] = data["basepair_counts"]["input"]
        d["bp_written"] = data["basepair_counts"]["output"]
        d["quality_trimmed"] = data["basepair_counts"]["quality_trimmed"]
        d["r_processed"] = data["read_counts"]["input"]
        d["r1_with_adapters"] = data["read_counts"]["read1_with_adapter"]
        d["r2_with_adapters"] = data["read_counts"]["read2_with_adapter"]
        d["r_too_short"] = data["read_counts"]["filtered"]["too_short"]
        d["r_too_long"] = data["read_counts"]["filtered"]["too_long"]
        d["r_too_many_N"] = data["read_counts"]["filtered"]["too_many_n"]
        d["r_written"] = data["read_counts"]["output"]
        d = {k: v for k, v in d.items() if v is not None}
        self.cutadapt_data[s_name] = d

        for end, end_key in [("5", "five_prime_end"), ("3", "three_prime_end")]:
            if end not in self.cutadapt_length_counts:
                self.cutadapt_length_counts[end] = dict()
                self.cutadapt_length_exp[end] = dict()
                self.cutadapt_length_obsexp[end] = dict()
            if s_name not in self.cutadapt_length_counts[end]:
                self.cutadapt_length_counts[end][s_name] = dict()
                self.cutadapt_length_exp[end][s_name] = dict()
                self.cutadapt_length_obsexp[end][s_name] = dict()

            for read in ["read1", "read2"]:
                if f"input_{read}" not in data["basepair_counts"]:
                    continue
                for adapter_data in data[f"adapters_{read}"]:
                    end_data = adapter_data.get(end_key)
                    if end_data:
                        for trimmed_length in end_data["trimmed_lengths"]:
                            length = trimmed_length["len"]
                            if length not in self.cutadapt_length_counts[end][s_name]:
                                self.cutadapt_length_counts[end][s_name][length] = 0
                                self.cutadapt_length_exp[end][s_name][length] = 0
                                self.cutadapt_length_obsexp[end][s_name][length] = 0
                            self.cutadapt_length_counts[end][s_name][length] += trimmed_length["counts"][0]
                            self.cutadapt_length_exp[end][s_name][length] += trimmed_length["expect"]
                            if trimmed_length["expect"] > 0:
                                self.cutadapt_length_obsexp[end][s_name][length] += (
                                    trimmed_length["counts"][0] / trimmed_length["expect"]
                                )

    def parse_log(self, f):
        """Go through log file looking for cutadapt output"""
        regexes = {
            "1.7": {
                "bp_processed": r"Total basepairs processed:\s*([\d,]+) bp",
                "bp_written": r"Total written \(filtered\):\s*([\d,]+) bp",
                "quality_trimmed": r"Quality-trimmed:\s*([\d,]+) bp",
                "r_processed": r"Total reads processed:\s*([\d,]+)",
                "pairs_processed": r"Total read pairs processed:\s*([\d,]+)",
                "r_with_adapters": r"Reads with adapters:\s*([\d,]+)",
                "r1_with_adapters": r"Read 1 with adapter:\s*([\d,]+)",
                "r2_with_adapters": r"Read 2 with adapter:\s*([\d,]+)",
                "r_too_short": r"Reads that were too short:\s*([\d,]+)",
                "pairs_too_short": r"Pairs that were too short:\s*([\d,]+)",
                "r_too_long": r"Reads that were too long:\s*([\d,]+)",
                "pairs_too_long": r"Pairs that were too long:\s*([\d,]+)",
                "r_too_many_N": r"Reads with too many N:\s*([\d,]+)",
                "pairs_too_many_N": r"Pairs with too many N:\s*([\d,]+)",
                "r_written": r"Reads written \(passing filters\):\s*([\d,]+)",
                "pairs_written": r"Pairs written \(passing filters\):\s*([\d,]+)",
            },
            "1.6": {
                "r_processed": r"Processed reads:\s*([\d,]+)",
                "bp_processed": r"Processed bases:\s*([\d,]+) bp",
                "r_trimmed": r"Trimmed reads:\s*([\d,]+)",
                "quality_trimmed": r"Quality-trimmed:\s*([\d,]+) bp",
                "bp_trimmed": r"Trimmed bases:\s*([\d,]+) bp",
                "too_short": r"Too short reads:\s*([\d,]+)",
                "too_long": r"Too long reads:\s*([\d,]+)",
            },
        }
        s_name = None
        end = "default"
        cutadapt_version = None
        parsing_version = "1.7"
        log_section = None
        path = os.path.join(f["root"], f["fn"])
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = iter(fh.readlines())
        for line in lines:
            # New log starting
            if "This is cutadapt" in line or "cutadapt version" in line:
                s_name = None
                end = "default"
                cutadapt_version = None
                c_version = re.match(r"This is cutadapt ([\d\.]+)", line)
                if c_version:
                    cutadapt_version = c_version.group(1)
                    try:
                        assert version.parse(c_version.group(1)) <= version.parse("1.6")
                        parsing_version = "1.6"
                    except Exception:
                        parsing_version = "1.7"
                c_version_old = re.match(r"cutadapt version ([\d\.]+)", line)
                if c_version_old:
                    cutadapt_version = c_version_old.group(1)
                    # The pattern "cutadapt version XX" is only pre-1.6
                    parsing_version = "1.6"
            # Get sample name from end of command line params
            cl_pref = "Command line parameters: "
            if line.startswith(cl_pref):
                input_fqs = []
                args = shlex.split(line[len(cl_pref) :])
                for i, x in enumerate(args):
                    if (
                        not x.startswith("-")
                        and x.endswith((".fastq", ".fq", ".gz", ".dat"))
                        and (i == 0 or args[i - 1] not in ["-o", "-p", "--output", "--paired-output"])
                    ):
                        input_fqs.append(x)
                if input_fqs:
                    s_name = self.clean_s_name(input_fqs, f)
                else:
                    # Manage case where sample name is '-' (reading from stdin)
                    s_name = f["s_name"]

                if s_name in self.cutadapt_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.cutadapt_data[s_name] = dict()
                if cutadapt_version:
                    self.cutadapt_data[s_name]["cutadapt_version"] = cutadapt_version

            if s_name is not None:
                # Add version info to module
                if cutadapt_version is not None:
                    self.add_software_version(cutadapt_version, s_name)

                self.add_data_source(f, s_name)

                # Search regexes for overview stats
                for k, r in regexes[parsing_version].items():
                    match = re.search(r, line)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(",", ""))

                # Starting a new section
                if "===" in line:
                    log_section = line.strip().strip("=").strip()

                # Detect whether 3' or 5'
                end_regex = re.search(r"Type: regular (\d)'", line)
                if end_regex:
                    end = end_regex.group(1)

                if "Overview of removed sequences" in line:
                    if "' end" in line:
                        res = re.search(r"(\d)' end", line)
                        end = res.group(1)

                    # Initialise dictionaries for length data if not already done
                    if end not in self.cutadapt_length_counts:
                        self.cutadapt_length_counts[end] = dict()
                        self.cutadapt_length_exp[end] = dict()
                        self.cutadapt_length_obsexp[end] = dict()

                # Histogram showing lengths trimmed
                if "length" in line and "count" in line and "expect" in line:
                    plot_sname = s_name
                    if log_section is not None:
                        plot_sname = f"{s_name} - {log_section}"
                    self.cutadapt_length_counts[end][plot_sname] = dict()
                    self.cutadapt_length_exp[end][plot_sname] = dict()
                    self.cutadapt_length_obsexp[end][plot_sname] = dict()

                    # Nested loop to read this section while the regex matches
                    for line2 in lines:
                        r_seqs = re.search(r"^(\d+)\s+(\d+)\s+([\d\.]+)", line2)
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

    def transform_trimming_length_data_for_plot(self):
        """Check if we parsed double ended data and transform it accordingly"""
        # Sanity check
        if (
            self.cutadapt_length_counts.keys() != self.cutadapt_length_exp.keys()
            or self.cutadapt_length_exp.keys() != self.cutadapt_length_obsexp.keys()
        ):
            log.error("Something went wrong...")
            log.debug("Keys in trimmed length data differed")
            raise ModuleNoSamplesFound

        if len(self.cutadapt_length_counts["default"]) == 0:
            self.cutadapt_length_counts.pop("default")
            self.cutadapt_length_exp.pop("default")
            self.cutadapt_length_obsexp.pop("default")

    def cutadapt_general_stats_table(self):
        """Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "percent_trimmed": {
                "title": "% BP Trimmed",
                "description": "% Total Base Pairs trimmed",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlBu-rev",
            }
        }
        self.general_stats_addcols(self.cutadapt_data, headers)

    def cutadapt_filtered_barplot(self):
        """Bar plot showing proportion of reads trimmed"""

        pconfig = {"id": "cutadapt_filtered_reads_plot", "title": "Cutadapt: Filtered Reads", "ylab": "Counts"}

        # We just use all categories. If a report is generated with a mixture
        # of SE and PE data then this means quite a lot of categories.
        # Usually, only a single data type is used though - in that case
        # any categories with 0 across all samples will be ignored.
        cats = {
            "pairs_written": {"name": "Pairs passing filters"},
            "r_written": {"name": "Reads passing filters"},
            "pairs_too_short": {"name": "Pairs that were too short"},
            "r_too_short": {"name": "Reads that were too short"},
            "pairs_too_long": {"name": "Pairs that were too long"},
            "r_too_long": {"name": "Reads that were too long"},
            "pairs_too_many_N": {"name": "Pairs with too many N"},
            "r_too_many_N": {"name": "Reads with too many N"},
            "pairs_filtered_unexplained": {"name": "Filtered pairs (uncategorised)"},
            "r_filtered_unexplained": {"name": "Filtered reads (uncategorised)"},
        }

        self.add_section(
            name="Filtered Reads",
            anchor="cutadapt_filtered_reads",
            description="This plot shows the number of reads (SE) / pairs (PE) removed by Cutadapt.",
            plot=bargraph.plot(self.cutadapt_data, cats, pconfig),
        )

    def cutadapt_length_trimmed_plot(self):
        """Generate the trimming length plot"""
        ends = [end for end, d_by_s in self.cutadapt_length_counts.items() if any(d_by_s.values())]
        ends = [x for x in ["default", "5", "3"] if x in ends]  # enforcing order
        for end in ends:
            pconfig = {
                "id": f"cutadapt_trimmed_sequences_plot_{end}",
                "title": "Cutadapt: Lengths of Trimmed Sequences{}".format(
                    "" if end == "default" else f" ({end}' end)"
                ),
                "ylab": "Counts",
                "xlab": "Length Trimmed (bp)",
                "ymin": 0,
                "tt_label": "<b>{point.x} bp trimmed</b>: {point.y:.0f}",
                "xsuffix": " bp",
                "data_labels": [
                    {"name": "Counts", "ylab": "Count"},
                    {"name": "Obs/Exp", "ylab": "Observed / Expected"},
                ],
            }

            self.add_section(
                name="Trimmed Sequence Lengths{}".format("" if end == "default" else f" ({end}')"),
                anchor=f"cutadapt_trimmed_sequences{'' if end == 'default' else f'_{end}'}",
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
