#!/usr/bin/env python

""" MultiQC module to parse output from Freyja """

import logging
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
            name="Freyja",
            anchor="freyja",
            href="https://github.com/andersen-lab/Freyja",
            info="Recover relative lineage abundances from mixed SARS-CoV-2 samples.",
            doi="10.1038/s41586-022-05049-6",
        )

        # To store the summary data
        self.freyja_data = dict()

        # Parse the output files
        self.parse_summ_files()

        # Remove filtered samples
        self.freyja_data = self.ignore_samples(self.freyja_data)

        # Let MultiQC know this module found no data
        if len(self.freyja_data) == 0:
            raise UserWarning

        log.info(f"Found {len(self.freyja_data)} reports")
        self.write_data_file(self.freyja_data, "multiqc_freyja")

        self.general_stats_cols()
        self.add_freyja_section()

    def parse_summ_files(self):
        """
        Parse the summary file.
        Freyja has multiple summary files, but we only need to parse the one from the demix command.
        More specifically, we only need the line that starts with "summarized".
        ...
        summarized	[('BQ.1*', 0.983), ('Omicron', 0.011), ('key', value)]
        ...
        """
        for f in self.find_log_files("freyja", filehandles=True):
            s_name = f["s_name"]
            # Read the statistics from file
            d = {}
            for line in f["f"]:
                try:
                    if line.startswith("summarized"):
                        summarized_line = line
                        summarized_line = summarized_line.strip().split("\t")[1]
                        d = eval(
                            summarized_line
                        )  # Make sure no input is corrupted and does not contain any malicious code
                        d = dict(d)
                except ValueError:
                    pass

            # Percentages don't always add up to 1, show a warning if this is the case
            if sum(d.values()) != 1:
                log.warning(f"Freyja {s_name}: percentages don't sum to 1")

            # There is no sample name in the log, so we use the root of the
            # file as sample name (since the filename is always stats.dat
            if s_name in self.freyja_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.freyja_data[s_name] = d
            self.add_data_source(f, s_name)

    def general_stats_cols(self):
        """Add a single column displaying the most abundant lineage to the General Statistics table"""
        top_lineage_dict = {}
        for s_name, sub_dict in self.freyja_data.items():
            top_lineage = max(sub_dict, key=sub_dict.get)
            top_lineage_value = sub_dict[top_lineage]
            top_lineage_dict[s_name] = {
                "Top_lineage_freyja": top_lineage,
                "Top_lineage_freyja_percentage": top_lineage_value,
            }

        headers = OrderedDict()
        headers["Top_lineage_freyja"] = {
            "title": "Top lineage (Freyja)",
            "description": "The most abundant lineage in the sample",
            "scale": "RdYlGn-rev",  # Not sure if this is the best scale
        }
        headers["Top_lineage_freyja_percentage"] = {
            "title": "Top lineage (Freyja) %",
            "description": "The percentage of the most abundant lineage in the sample",
            "max": 100,
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * 100,
            "suffix": "%",
        }

        self.general_stats_addcols(top_lineage_dict, headers)

    def add_freyja_section(self):
        cats = list()
        pconfig = {
            "id": "Freyja_plot",
            "title": "Freyja: Top lineages",
            "ylab": "relative abundance",
            "yCeiling": 1,
            "cpswitch": False,
            "cpswitch_c_active": False,
        }

        rank_cats = OrderedDict()

        rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
        cats.append(rank_cats)

        self.add_section(
            name="Freyja Summary",
            anchor="freyja-summary",
            description="""
                Relative lineage abundances from mixed SARS-CoV-2 samples. Hover over the column headers for descriptions and click _Help_ for more in-depth documentation. 
                """,
            helptext="""
                The graph denotes a sum of all lineage abundances in a particular WHO designation , otherwise they are grouped into "Other".
                Lineages abundances are calculated as the number of reads that are assigned to a particular lineage. 
                Lineages and their corresponding abundances are summarized by constellation. 

                > **Note**: Lineage designation is based on the used WHO nomenclature, which could vary over time. 
                """,
            plot=bargraph.plot(self.freyja_data, None, pconfig),
        )
