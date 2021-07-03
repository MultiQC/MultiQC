#!/usr/bin/env python

""" MultiQC module to parse output from NanoStat """

from __future__ import print_function
from collections import OrderedDict
import logging
import jinja2

from multiqc import config
from multiqc.utils import mqc_colour
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """NanoStat module"""

    _KEYS_NUM = [
        "Active channels",
        "Number of reads",
        "Total bases",
        "Total bases aligned",
        "Read length N50",
        "Mean read length",
        "Median read length",
        "Median read quality",
        "Mean read quality",
        "Average percent identity",
        "Median percent identity",
    ]

    _KEYS_READ_Q = [
        "&gt;Q5",
        "&gt;Q7",
        "&gt;Q10",
        "&gt;Q12",
        "&gt;Q15",
    ]
    _stat_types = ("aligned", "seq summary", "unrecognized")

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="NanoStat",
            anchor="nanostat",
            href="https://github.com/wdecoster/nanostat/",
            info="various statistics from a long read sequencing dataset in fastq, bam or sequencing summary format.",
        )

        # Find and load any NanoStat reports
        self.nanostat_data = dict()
        self.has_aligned = False
        self.has_seq_summary = False
        for f in self.find_log_files("nanostat", filehandles=True):
            self.parse_nanostat_log(f)

        # Filter to strip out ignored sample names
        self.nanostat_data = self.ignore_samples(self.nanostat_data)

        if len(self.nanostat_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.nanostat_data)))

        # Write parsed report data to a file
        self.write_data_file(self.nanostat_data, "multiqc_nanostat")

        # Stats Tables
        if self.has_aligned:
            self.nanostat_stats_table("aligned")
        if self.has_seq_summary:
            self.nanostat_stats_table("seq summary")

        # Quality distribution Plot
        self.reads_by_quality_plot()

    def parse_nanostat_log(self, f):
        """Parse output from NanoStat

        Note: Tool can be run in two different modes, giving two variants to the output.
        To avoid overwriting keys from different modes, keys are given a suffix.
        """

        nano_stats = {}
        for line in f["f"]:

            line = jinja2.escape(line)
            parts = line.strip().split(":")
            if len(parts) == 0:
                continue

            key = parts[0]

            if key in self._KEYS_NUM:
                val = float(parts[1].replace(",", ""))
                nano_stats[key] = val
            elif key in self._KEYS_READ_Q:
                # Number of reads above Q score cutoff
                val = int(parts[1].strip().split()[0])
                nano_stats[key] = val

        if "Total bases aligned" in nano_stats:
            stat_type = "aligned"
            self.has_aligned = True
        elif "Active channels" in nano_stats:
            stat_type = "seq summary"
            self.has_seq_summary = True
        else:
            log.debug(f"Did not recognise NanoStat file '{f['fn']}' - skipping")
            return

        out_d = {f"{k}_{stat_type}": v for k, v in nano_stats.items()}

        # Warn if we find overlapping data for the same sample
        if f["s_name"] in self.nanostat_data:
            # Only if the same has some keys in common
            if not set(self.nanostat_data[f["s_name"]].keys()).isdisjoint(out_d.keys()):
                log.debug("Duplicate sample data found! Overwriting: {}".format(f["s_name"]))

        self.nanostat_data.setdefault(f["s_name"], {}).update(out_d)

        self.add_data_source(f)

    def nanostat_stats_table(self, stat_type):
        """Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report"""

        headers_base = OrderedDict()
        headers_base["Active channels"] = {
            "title": "Active channels",
            "description": "Active channels",
            "scale": "Greens",
            "format": "{:,.0f}",
        }
        headers_base["Median read length"] = {
            "title": f"Median length",
            "description": f"Median read length (bp)",
            "suffix": " bp",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
            "scale": "BuPu",
        }
        headers_base["Mean read length"] = {
            "title": f"Mean length",
            "description": f"Mean read length (bp)",
            "suffix": " bp",
            "scale": "Purples",
            "format": "{:,.0f}",
            "shared_key": "nucleotides",
            "hidden": True,
        }
        headers_base["Read length N50"] = {
            "title": "Read N50",
            "description": "Read length N50",
            "format": "{:,.0f}",
            "suffix": " bp",
            "scale": "RdPu",
        }
        headers_base["Median read quality"] = {
            "title": "Median Qual",
            "description": "Median read quality (Phred scale)",
            "shared_key": "phred_score",
            "scale": "RdYlGn",
        }
        headers_base["Mean read quality"] = {
            "title": "Mean Qual",
            "description": "Mean read quality (Phred scale)",
            "scale": "PiYG",
            "shared_key": "phred_score",
            "hidden": True,
        }
        headers_base["Median percent identity"] = {
            "title": "Median Identity",
            "description": "Median percent identity",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlBu",
            "shared_key": "percent_identity",
        }
        headers_base["Average percent identity"] = {
            "title": "Mean Identity",
            "description": "Average percent identity",
            "max": 100,
            "suffix": "%",
            "scale": "Spectral",
            "shared_key": "percent_identity",
            "hidden": True,
        }
        headers_base["Number of reads"] = {
            "title": f"# Reads ({config.long_read_count_prefix})",
            "description": f"Number of reads ({config.long_read_count_desc})",
            "modify": lambda x: x * config.long_read_count_multiplier,
            "shared_key": "long_read_count",
            "scale": "YlGn",
        }
        headers_base["Total bases"] = {
            "title": f"Total Bases ({config.base_count_prefix})",
            "description": f"Total bases ({config.base_count_desc})",
            "modify": lambda x: x * config.base_count_multiplier,
            "shared_key": "base_count",
            "scale": "BrBG",
        }
        headers_base["Total bases aligned"] = {
            "title": f"Aligned Bases ({config.base_count_prefix})",
            "description": f"Total bases aligned ({config.base_count_desc})",
            "modify": lambda x: x * config.base_count_multiplier,
            "shared_key": "base_count",
            "scale": "PuOr",
        }

        # Add the stat_type suffix
        headers = OrderedDict()
        for k in headers_base:
            key = f"{k}_{stat_type}"
            headers[key] = headers_base.get(k, dict()).copy()

        # Table config
        table_config = {
            "namespace": "NanoStat",
            "id": "nanostat_{}_stats_table".format(stat_type.replace(" ", "_")),
            "table_title": f"NanoStat {stat_type}",
        }

        # Add the report section
        description = ""
        if stat_type == "aligned":
            description = "NanoStat statistics from FastQ, FASTA or BAM files."
        if stat_type == "seq summary":
            description = "NanoStat statistics from albacore or guppy summary files."
        self.add_section(
            name="{} stats".format(stat_type.replace("_", " ").capitalize()),
            anchor="nanostat_{}_stats".format(stat_type.replace(" ", "_")),
            description=description,
            plot=table.plot(self.nanostat_data, headers, table_config),
        )

    def reads_by_quality_plot(self):
        """Make the HighCharts HTML to plot the reads by quality"""

        def _get_total_reads(data_dict):
            stat_type = self._stat_types[0]
            for stat_type in self._stat_types:
                total_key = f"Number of reads_{stat_type}"
                if total_key in data_dict:
                    return data_dict[total_key], stat_type
            return None, None

        bar_data = {}
        stat_type = "unrecognized"
        # Order of keys, from >Q5 to >Q15
        _range_names = {
            "&gt;Q5": "&lt;Q5",
            "&gt;Q7": "Q5-7",
            "&gt;Q10": "Q7-10",
            "&gt;Q12": "Q10-12",
            "&gt;Q15": "Q12-15",
            "rest": "&gt;Q15",
        }
        for s_name, data_dict in self.nanostat_data.items():
            reads_total, stat_type = _get_total_reads(data_dict)
            if s_name in bar_data and stat_type == "aligned":
                log.debug("Sample '{s_name}' duplicated in the quality plot - ignoring aligned data")
                continue
            elif s_name in bar_data and stat_type == "seq summary":
                log.debug("Sample '{s_name}' duplicated in the quality plot - overwriting with seq summary data")
            bar_data[s_name] = {}

            prev_reads = reads_total
            for k, range_name in _range_names.items():
                if k != "rest":
                    data_key = f"{k}_{stat_type}"
                    reads_gt = data_dict[data_key]

                    bar_data[s_name][range_name] = prev_reads - reads_gt

                    if bar_data[s_name][range_name] < 0:
                        log.error(f"Error on {s_name} {range_name} {data_key} . Negative number of reads")
                    prev_reads = reads_gt
                else:
                    data_key = f"&gt;Q15_{stat_type}"
                    bar_data[s_name][range_name] = data_dict[data_key]

        cats = OrderedDict()
        keys = reversed(list(_range_names.values()))
        colours = mqc_colour.mqc_colour_scale("RdYlGn-rev", 0, len(_range_names))
        for idx, k in enumerate(keys):
            cats[k] = {"name": "Reads " + k, "color": colours.get_colour(idx, lighten=1)}

        # Config for the plot
        config = {
            "id": "nanostat_quality_dist",
            "title": "NanoStat: Reads by quality",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        # Add the report section
        self.add_section(
            name="Reads by quality",
            anchor=f"nanostat_read_qualities",
            description="Read counts categorised by read quality (phred score).",
            helptext="""
                Sequencing machines assign each generated read a quality score using the
                [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score).
                The phred score represents the liklelyhood that a given read contains errors.
                So, high quality reads have a high score.

                Data may come from NanoPlot reports generated with sequencing summary files or alignment stats.
                If a sample has data from both, the sequencing summary is preferred.
            """,
            plot=bargraph.plot(bar_data, cats, config),
        )
