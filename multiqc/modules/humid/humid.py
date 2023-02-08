""" MultiQC module to parse output from Lima """

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="HUMID",
            anchor="humid",
            href="https://github.com/jfjlaros/HUMID",
            info=" is a fast, reference free tool to remove (UMI) duplicates from sequencing data",
            # No publication / DOI // doi=
        )

        # To store the summary data
        self.humid = dict()

        # Parse the output files
        self.parse_stat_files()

        # Remove filtered samples
        self.humid = self.ignore_samples(self.humid)

        # Let MultiQC know this module found no data
        if not self.humid:
            raise UserWarning

        log.info(f"Found {len(self.humid)} reports")
        self.write_data_file(self.humid, "multiqc_humid")

        self.add_general_stats()
        self.add_humid_section()

    def parse_stat_files(self):
        for f in self.find_log_files("humid", filehandles=True):
            # There is no sample name in the log, so we use the root of the
            # file as sample name (since the filename is always stats.dat
            s_name = self.clean_s_name(f["root"], f)

            # Read the statistics from file
            d = {}
            for line in f["f"]:
                try:
                    field, value = line.strip().split(": ")
                    d[field] = int(value)
                except ValueError:
                    pass

            try:
                # Calculate additional statistics
                d["filtered"] = d["total"] - d["usable"]
                d["duplicates"] = d["total"] - d["clusters"] - d["filtered"]

                # Make sure we only return data that makes sense
                if not sum(d[field] for field in ["duplicates", "clusters", "filtered"]) == d["total"]:
                    log.warning(f"HUMID stats looked wrong, skipping: {s_name}")
                    continue
            except KeyError as e:
                log.warning(f"Expected HUMID keys missing {e}, skipping sample: {s_name}")
                continue

            # Got this far, data must be good
            if s_name in self.humid:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.humid[s_name] = d
            self.add_data_source(f, s_name)

    def add_general_stats(self):
        # Add the number of unique reads (=clusters) to the general statistics
        # report
        data = {k: {"uniq": v["clusters"]} for k, v in self.humid.items()}
        headers = OrderedDict()
        headers["uniq"] = {
            "title": f"{config.read_count_prefix} Unique",
            "description": f"Number of unique reads after UMI deduplication ({config.read_count_desc})",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        self.general_stats_addcols(data, headers)

    def add_humid_section(self):
        # The values we want to plot (add to the total number of reads)
        cats = OrderedDict()
        cats["clusters"] = {"name": "Unique reads"}
        cats["duplicates"] = {"name": "Duplicate reads"}
        cats["filtered"] = {"name": "Filtered reads"}

        # Bargraph configuration
        config = {
            "id": "humid-bargraph",
            "title": "HUMID: Deduplication results",
            "ylab": "Number of reads",
            "hide_zero_cats": False,
        }
        self.add_section(
            name="Duplication Summary",
            anchor="humid-section",
            description="""
                Duplication statistics per sample. Every read in the
                input data has been assigned to one of the three categories
                shown here.
                """,
            helptext="""
                - **Unique reads** are the reads that are left over after deduplication.
                - **Duplicate reads** are the reads that were determined to be duplicates of the **Unique reads**.
                - **Filtered reads** were reads that could not be analysed, due to N nucleotides, or because the UMI or sequences were too short to use.
                """,
            plot=bargraph.plot(self.humid, cats, config),
        )
