#!/usr/bin/env python

""" MultiQC module to parse output from NanoStat """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ NanoStat module """

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
    _key_fmt = "{} ({})"

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
        for f in self.find_log_files("nanostat", filehandles=True):
            self.parse_nanostat_log(f)

        # Filter to strip out ignored sample names
        self.nanostat_data = self.ignore_samples(self.nanostat_data)

        if len(self.nanostat_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.nanostat_data)))

        # Write parsed report data to a file
        self.write_data_file(self.nanostat_data, "multiqc_nanostat")

        # Basic Stats Table
        self.nanostat_general_stats_table()

        # Quality distribution Plot
        self.add_section(plot=self.nanostat_alignment_plot())

    def parse_nanostat_log(self, logfile):
        import jinja2

        s_name = logfile["s_name"]
        nano_stats = {}
        for line in logfile["f"]:

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
        elif "Active channels" in nano_stats:
            stat_type = "seq summary"
        else:
            stat_type = "unrecognized"
            log.debug("%s datatype with key %s", stat_type, str(logfile))
            log.debug(str(nano_stats))

        out_d = {self._key_fmt.format(k, stat_type): v for k, v in nano_stats.items()}

        # Warn if we find overlapping data for the same sample
        if logfile["s_name"] in self.nanostat_data and not set(
            self.nanostat_data[s_name].keys()
        ).isdisjoint(out_d.keys()):
            log.debug(
                "Duplicate sample data found! Overwriting: {}".format(logfile["s_name"])
            )

        self.nanostat_data.setdefault(s_name, dict()).update(out_d)

        self.add_data_source(logfile, s_name)

    def nanostat_general_stats_table(self):
        """Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report"""

        headers_base = {
            "Average percent identity": {
                "title": "Mean Identity",
                "description": "Average percent identity",
                "max": 100,
                # "min": 0,
                "suffix": "%",
                "scale": "YlGn",
                "hidden": True,
            },
            "Active channels": {
                "title": "Active channels",
                "description": "Active channels",
                "format": "{:,g}",
                # "min": 0,
                "scale": "YlGn",
            },
            "Mean read length": {
                "title": "Mean length",
                "description": "Mean read length",
                # "min": 0,
                "format": "{:,g}",
                "suffix": "bp",
                "scale": "YlGn",
                "hidden": True,
            },
            "Mean read quality": {
                "title": "Mean Qual",
                "description": "Mean read quality (Phred scale)",
                # "min": 0,
                "scale": "YlGn",
                "hidden": True,
            },
            "Median percent identity": {
                "title": "Median Identity",
                "description": "Median percent identity",
                "max": 100,
                # "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            },
            "Median read length": {
                "title": "Median length",
                "description": "Median read length",
                # "min": 0,
                "suffix": "bp",
                "format": "{:,g}",
                "scale": "YlGn",
                "hidden": True,
            },
            "Median read quality": {
                "title": "Median Qual",
                "description": "Median read quality (Phred scale)",
                # "min": 0,
                "scale": "YlGn",
            },
            "Number of reads": {
                "title": "# reads",
                "description": "Number of reads",
                "format": "{:,g}",
                # "min": 0,
                "scale": "YlGn",
            },
            "Read length N50": {
                "title": "Read N50",
                "description": "Read length N50",
                "format": "{:,g}",
                # "min": 0,
                "suffix": "bp",
                "scale": "YlGn",
            },
            "Total bases": {
                "title": "Total Bases",
                "description": "Total bases",
                "format": "{:,g}",
                # "min": 0,
                "suffix": "bp",
                "scale": "YlGn",
            },
            "Total bases aligned": {
                "title": "Aligned Bases",
                "description": "Total bases aligned",
                "format": "{:,g}",
                # "min": 0,
                "suffix": "bp",
                "scale": "YlGn",
            },
        }
        # Find all columns to add to the table"
        header_keys = set()
        for d in self.nanostat_data.values():
            header_keys.update({k for k in d.keys()})

        # table_keys = list(next(iter(self.nanostat_data.values())).keys())
        default_stat_source = self._stat_types[0]

        headers = dict()
        for placement_idx, k in enumerate(self._KEYS_NUM + self._KEYS_READ_Q):
            first_type_key = self._key_fmt.format(k, default_stat_source)
            for stat_type_idx, stat_type in enumerate(self._stat_types):
                key = self._key_fmt.format(k, stat_type)
                if key in header_keys:
                    headers[key] = headers_base.get(k, dict()).copy()
                    headers[key]["placement"] = 1000.0 + placement_idx

                    if "title" in headers[key]:
                        headers[key]["title"] += f" ({stat_type})"
                    # Hide columns from sequencing_summary if alignment summary available.
                    visible = (stat_type_idx == 0) or (first_type_key not in headers)
                    visible &= not key.startswith("&gt;")
                    headers[key]["hidden"] = headers[key].get("hidden", False) or (
                        visible == False
                    )

        self.general_stats_addcols(self.nanostat_data, headers)

    def nanostat_alignment_plot(self):
        """ Make the HighCharts HTML to plot the alignment rates """

        def _get_total_reads(data_dict):
            stat_type = self._stat_types[0]
            for stat_type in self._stat_types:
                total_key = self._key_fmt.format("Number of reads", stat_type)
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
            bar_data[s_name] = {}
            reads_total, stat_type = _get_total_reads(data_dict)

            prev_reads = reads_total
            for k, range_name in _range_names.items():
                if k != "rest":
                    data_key = self._key_fmt.format(k, stat_type)
                    reads_gt = data_dict[data_key]

                    bar_data[s_name][range_name] = prev_reads - reads_gt

                    assert (
                        bar_data[s_name][range_name] >= 0
                    ), f"Error on {s_name} {range_name} {data_key} . Negative number of reads"
                    prev_reads = reads_gt
                else:
                    data_key = self._key_fmt.format("&gt;Q15", stat_type)
                    bar_data[s_name][range_name] = data_dict[data_key]

        keys = OrderedDict()

        for k in _range_names.values():
            keys[k] = {"name": "Reads " + k}

        # Config for the plot
        config = {
            "id": "nanostat_quality_dist",
            "title": f"NanoStat: Reads by quality ({stat_type})",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(bar_data, keys, config)
