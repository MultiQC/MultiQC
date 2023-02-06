""" MultiQC module to parse log output from Xenome Classify """


import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Xenome",
            anchor="xenome",
            href="https://github.com/data61/gossamer/blob/master/docs/xenome.md",
            info="is a tool for classifying reads from xenograft sources.",
            doi="doi.org/10.1093/bioinformatics/bts236",
        )

        # Parse logs
        # Find and load any Bowtie reports
        self.xenome_data = dict()
        for f in self.find_log_files("xenome"):
            self.parse_xenome_logs(f)

        # Filter to strip out ignored sample names
        self.xenome_data = self.ignore_samples(self.xenome_data)

        if len(self.xenome_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.xenome_data)))

        # Write parsed report data to a file
        self.write_data_file(self.xenome_data, "multiqc_xenome")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.xenome_general_stats_table()

        # Alignment Rate Plot
        self.xenome_stats_plot()

    def parse_xenome_logs(self, f):
        s_name = f["s_name"].replace("_xenome_stats", "")
        parsed_data = {}
        regexes = {
            "human_reads": r"(\d+)\s+\d+.\d+\s+human",
            "mouse_reads": r"(\d+)\s+\d+.\d+\s+mouse",
            "both_human_mouse_reads": r"(\d+)\s+\d+.\d+\s+both",
            "neither_human_mouse_reads": r"(\d+)\s+\d+.\d+\s+neither",
            "ambiguous_reads": r"(\d+)\s+\d+.\d+\s+ambiguous",
        }

        for l in f["f"].splitlines():
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = float(match.group(1))

        if len(parsed_data) > 0:
            if s_name in self.xenome_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.xenome_data[s_name] = parsed_data

    def xenome_general_stats_table(self):
        """Take the parsed stats from the Xenome log and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["human_reads"] = {
            "title": "Number Human Reads",
            "description": "The number of human reads in the sample",
            "min": 0,
            "scale": False,
        }
        headers["mouse_reads"] = {
            "title": "Number Mouse Reads",
            "description": "The number of human reads in the sample",
            "min": 0,
            "scale": False,
        }
        self.general_stats_addcols(self.xenome_data, headers)

    def xenome_stats_plot(self):
        """Make the HighCharts HTML to plot the alignment rates"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["human_reads"] = {"color": "#000034", "name": "Human"}
        keys["mouse_reads"] = {"color": "#328AE2", "name": "Mouse"}
        keys["both_human_mouse_reads"] = {"color": "#7F9DA7", "name": "Both Species"}
        keys["neither_human_mouse_reads"] = {"color": "#C1B8AA", "name": "Neither Species"}
        keys["ambiguous_reads"] = {"color": "#DDD4C6", "name": "Ambiguous"}

        # Config for the plot
        config = {
            "id": "xenome_stats",
            "title": "Xenome: Classification Counts",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(
            description="This plot shows the number of reads classified by Xenome to human and/or mouse.",
            helptext="""
            There are 5 possible categories:
            * **Human**: Read was found only in human.
            * **Mouse**: Read was found only in mouse.
            * **Both Species**: Read was found in either mouse or human.
            * **Neither Species**: Read was found in neither mouse or human.
            * **Ambiguous**: Read origin could not be adequately determined.
            """,
            plot=bargraph.plot(self.xenome_data, keys, config),
        )
