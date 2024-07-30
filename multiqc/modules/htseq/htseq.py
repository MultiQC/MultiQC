import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    HTSeq is a general purpose Python package that provides infrastructure to
    process data from high-throughput sequencing assays. `htseq-count` is a tool
    that is part of the main HTSeq package - it takes a file with aligned sequencing
    reads, plus a list of genomic features and counts how many reads map to each feature.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HTSeq Count",
            anchor="htseq",
            target="HTSeq Count",
            href="https://htseq.readthedocs.io/en/master/htseqcount.html",
            info="Part of the HTSeq package: counts reads covering specified genomic features",
            doi="10.1093/bioinformatics/btu638",
        )

        # Find and load any HTSeq Count reports
        self.htseq_data = dict()
        self.htseq_keys = list()
        for f in self.find_log_files("htseq", filehandles=True):
            parsed_data = self.parse_htseq_report(f)
            if parsed_data is not None:
                self.htseq_data[f["s_name"]] = parsed_data
                self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.htseq_data = self.ignore_samples(self.htseq_data)

        if len(self.htseq_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.htseq_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.htseq_data, "multiqc_htseq")

        # Basic Stats Table
        self.htseq_stats_table()

        # Assignment bar plot
        self.add_section(plot=self.htseq_counts_chart())

    @staticmethod
    def parse_htseq_report(f):
        """Parse the HTSeq Count log file."""
        keys = ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
        parsed_data = dict()
        assigned_counts = 0

        # HtSeq search pattern is just two tab-separated columns, which is not very specific.
        # Need to wrap in try-catch to catch potential weird files like parquet being matched.
        try:
            for line in f["f"]:
                s = line.split("\t")
                if s[0] in keys:
                    parsed_data[s[0][2:]] = int(s[-1])
                else:
                    try:
                        assigned_counts += int(s[-1])
                    except (ValueError, IndexError):
                        pass
        except UnicodeDecodeError:
            log.debug(f"Could not parse potential HTSeq Count file {f['fn']}")
            return None

        if len(parsed_data) > 0:
            parsed_data["assigned"] = assigned_counts
            parsed_data["total_count"] = sum([v for v in parsed_data.values()])
            try:
                parsed_data["percent_assigned"] = (
                    float(parsed_data["assigned"]) / float(parsed_data["total_count"])
                ) * 100.0
            except ZeroDivisionError:
                parsed_data["percent_assigned"] = 0
            return parsed_data
        return None

    def htseq_stats_table(self):
        """Take the parsed stats from the HTSeq Count report and add them to the
        basic stats table at the top of the report"""

        headers = {
            "percent_assigned": {
                "title": "% Assigned",
                "description": "% Assigned reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "assigned": {
                "title": f"{config.read_count_prefix} Assigned",
                "description": f"Assigned Reads ({config.read_count_desc})",
                "min": 0,
                "scale": "PuBu",
                "modify": lambda x: float(x) * config.read_count_multiplier,
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.htseq_data, headers)

    def htseq_counts_chart(self):
        """Make the HTSeq Count assignment rates plot"""
        cats = {
            "assigned": {"name": "Assigned"},
            "ambiguous": {"name": "Ambiguous"},
            "alignment_not_unique": {"name": "Alignment Not Unique"},
            "no_feature": {"name": "No Feature"},
            "too_low_aQual": {"name": "Too Low aQual"},
            "not_aligned": {"name": "Not Aligned"},
        }
        config = {
            "id": "htseq_assignment_plot",
            "title": "HTSeq: Count Assignments",
            "ylab": "# Reads",
            "hide_empty": False,
            "cpswitch_counts_label": "Number of Reads",
        }
        return bargraph.plot(self.htseq_data, cats, config)
