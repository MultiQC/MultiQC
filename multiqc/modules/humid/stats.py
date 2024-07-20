import logging
from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


def parse_reports(self):
    # To store the summary data
    self.stats = dict()

    # Parse the output files
    parse_stat_files(self)

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Remove filtered samples
    self.stats = self.ignore_samples(self.stats)

    if self.stats:
        self.write_data_file(self.stats, "multiqc_humid_stats")
        add_general_stats(self)
        add_stats_section(self)


def parse_stat_files(self):
    for f in self.find_log_files("humid/stats", filehandles=True):
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
        if s_name in self.stats:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.stats[s_name] = d
        self.add_data_source(f, s_name)


def add_general_stats(self):
    def add_unique_reads(self):
        """
        Add the number of reads after deduplication

        This corresponds to the number of clusters in HUMID
        """
        data = {k: {"uniq": v["clusters"]} for k, v in self.stats.items()}
        headers = dict()
        headers["uniq"] = {
            "title": f"{config.read_count_prefix} Unique Reads",
            "description": f"Reads remaining after deduplication ({config.read_count_desc})",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Blues",
        }
        self.general_stats_addcols(data, headers)

    def add_passing_deduplication(self):
        """
        Add the percentage of reads remaining after deduplication
        """
        data = {k: {"perc": (v["clusters"] / v["total"]) * 100} for k, v in self.stats.items()}
        headers = dict()
        headers["perc"] = {
            "title": "% Pass Dedup",
            "description": "% processed reads that passed deduplication",
            "suffix": "%",
            "min": 0,
            "max": 100,
            "scale": "RdBu",
        }
        self.general_stats_addcols(data, headers)

    add_unique_reads(self)
    add_passing_deduplication(self)


def add_stats_section(self):
    # The values we want to plot (add to the total number of reads)
    cats = dict()
    cats["clusters"] = {"name": "Unique reads"}
    cats["duplicates"] = {"name": "Duplicate reads"}
    cats["filtered"] = {"name": "Filtered reads"}

    # Bargraph configuration
    config = {
        "id": "humid-bargraph",
        "title": "HUMID: Deduplication results",
        "ylab": "Number of reads",
        "hide_empty": False,
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
        plot=bargraph.plot(self.stats, cats, config),
    )
