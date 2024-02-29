""" MultiQC module to parse logs from Skewer """


import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"skewer v([\d\.]+) \[.+\]"


class MultiqcModule(BaseMultiqcModule):
    """Skewer"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Skewer",
            anchor="skewer",
            href="https://github.com/relipmoc/skewer",
            info="is an adapter trimming tool specially designed for processing next-generation sequencing (NGS) paired-end sequences.",
            doi="10.1186/1471-2105-15-182",
        )

        self.skewer_data = dict()
        self.skewer_readlen_dist = dict()

        for f in self.find_log_files("skewer", filehandles=True):
            self.parse_skewer_log(f)

        # Filter to strip out ignored sample names
        self.skewer_data = self.ignore_samples(self.skewer_data)

        if len(self.skewer_data) == 0:
            raise ModuleNoSamplesFound

        headers = {
            "pct_trimmed": {
                "title": "% Trimmed",
                "description": "% of reads trimmed",
                "scale": "RdYlGn-rev",
                "max": 100,
                "min": 0,
                "suffix": "%",
            }
        }

        self.general_stats_addcols(self.skewer_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.skewer_data, "multiqc_skewer")
        self.write_data_file(self.skewer_readlen_dist, "multiqc_skewer_readlen_dist")

        # set the value 0 for every x where a given sample doens't have a value
        all_x_values = []
        for s_name in self.skewer_readlen_dist:
            for xval in self.skewer_readlen_dist[s_name]:
                all_x_values.append(xval)

        for s_name in self.skewer_readlen_dist:
            for xval in all_x_values:
                if xval not in self.skewer_readlen_dist[s_name]:
                    self.skewer_readlen_dist[s_name][xval] = 0.0

        # add the histogram to the report
        self.add_readlen_dist_plot()

        log.info(f"Found {len(self.skewer_data)} reports")

    def add_readlen_dist_plot(self):
        """Generate plot HTML for read length distribution plot."""
        pconfig = {
            "id": "skewer_read_length_histogram",
            "title": "Skewer: Read Length Distribution after trimming",
            "xDecimals": False,
            "ylab": "% of Reads",
            "xlab": "Read Length",
            "xmin": 0,
            "ymin": 0,
            "tt_label": "<b>{point.x}</b>: {point.y:.1f}%",
        }
        self.add_section(plot=linegraph.plot(self.skewer_readlen_dist, pconfig))

    def parse_skewer_log(self, f):
        """Go through log file looking for skewer output"""
        fh = f["f"]
        regexes = {
            "fq1": r"Input file:\s+(.+)",
            "fq2": r"Paired file:\s+(.+)",
            "r_processed": r"(\d+) read|reads pairs? processed",
            "r_short_filtered": r"(\d+) \(\s*\d+.\d+%\) short read",
            "r_empty_filtered": r"(\d+) \(\s*\d+.\d+%\) empty read",
            "r_avail": r"(\d+) \(\s*\d+.\d+%\) read",
            "r_trimmed": r"(\d+) \(\s*\d+.\d+%\) trimmed read",
            "r_untrimmed": r"(\d+) \(\s*\d+.\d+%\) untrimmed read",
        }
        regex_hist = r"\s?(\d+)\s+(\d+)\s+(\d+.\d+)%"

        data = dict()
        for k, v in regexes.items():
            data[k] = 0
        data["fq1"] = None
        data["fq2"] = None
        readlen_dist = dict()

        for line in fh:
            if line.startswith("skewer"):
                match = re.search(VERSION_REGEX, line)
                if match:
                    data["version"] = match.group(1)

            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    data[k] = match.group(1).replace(",", "")

            match = re.search(regex_hist, line)
            if match:
                read_length = int(match.group(1))
                pct_at_rl = float(match.group(3))
                readlen_dist[read_length] = pct_at_rl

        if data["fq1"] is not None:
            s_name = self.clean_s_name(data["fq1"], f)
            if s_name in self.skewer_readlen_dist:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.add_skewer_data(s_name, data, f)
            self.skewer_readlen_dist[s_name] = readlen_dist
            self.add_software_version(data.get("version"), s_name)

        if data["fq2"] is not None:
            s_name = self.clean_s_name(data["fq1"], f)
            if s_name in self.skewer_readlen_dist:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.add_skewer_data(s_name, data, f)
            self.skewer_readlen_dist[s_name] = readlen_dist
            self.add_software_version(data.get("version"), s_name)

    def add_skewer_data(self, s_name, data, f):
        stats = ["r_processed", "r_short_filtered", "r_empty_filtered", "r_avail", "r_trimmed", "r_untrimmed"]
        if s_name in self.skewer_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.skewer_data[s_name] = {}
        self.add_data_source(f, s_name)
        for k in stats:
            self.skewer_data[s_name][k] = int(data[k])

        self.skewer_data[s_name]["pct_avail"] = (
            100.0 * float(data["r_avail"]) / float(data["r_processed"]) if float(data["r_processed"]) else None
        )
        self.skewer_data[s_name]["pct_trimmed"] = (
            100.0 * float(data["r_trimmed"]) / float(data["r_avail"]) if float(data["r_avail"]) else None
        )
        self.skewer_data[s_name]["pct_untrimmed"] = (
            100.0 * float(data["r_untrimmed"]) / float(data["r_avail"]) if float(data["r_avail"]) else None
        )
