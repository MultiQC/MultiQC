import logging
import json

from typing import List
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics from the JSON log file (--summary or -s options). It parses and reports
    the number of input, output and removed sequences displays them in the General Stats table as well as a bar plot.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Deacon",
            anchor="deacon",
            href="https://github.com/bede/deacon",
            doi="https://doi.org/10.1101/2025.06.09.658732",
            info="Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching.",
        )

        self.deacon_data = dict()
        # self. makes deacon.data an instance variable, not a local for __init__

        # iterate through all detected Deacon log files
        # includes an open file handle in the returned dict under key "f"
        for f in self.find_log_files("deacon", filehandles=True):  # filehandles: opens the files
            try:
                # load JSON data from filehandles, if `f["f"]` is a file-like object
                data = json.load(f["f"])
            except Exception as e:  # catch parsing-errors in e
                log.warning(f"failed to parse data from {f['fn']} : {e}")  # {e} -> error message
                continue

            # check, if there is "version" in JSON report and contains "version"
            version = data.get("version", "")
            if not version.startswith("deacon"):  # if not, give out a warning
                log.warning(f"Skipping {f['fn']} : no deacon report or wrong version")
                continue

            self.add_data_source(f)

            sample_name = self.clean_s_name(
                f["s_name"], f
            )  # f["s_name"] -> samplename; self.clean_s_name -> normalizes the name to avoid duplicates

            # extract metrics for each sample
            self.deacon_data[sample_name] = {
                # "k" : data.get("k"), #k-mere length
                # "w" : data.get("w"), #wondow size
                # "abs_threshold" : data.get("abs_threshold"), #number of k-mers that count as hit
                # "rel_threshold" : data.get("rel_threshold"),
                # "prefix_length" : data.get("prefix_length"), #prefix length used for hashing/filtering
                "deplete": ("True" if data.get("deplete") else "False"),  # boolean, deplete mode (True/False)
                # "rename" : data.get("rename"), #boolean, renamed seqs after export (True/False)
                "seqs_in": data.get("seqs_in"),  # number of input seqs
                "seqs_out": data.get("seqs_out"),  # number of output seqs
                "seqs_out_proportion": data.get("seqs_out_proportion"),
                "seqs_removed": data.get("seqs_removed"),  # number of removed seqs
                "seqs_removed_proportion": data.get("seqs_removed_proportion"),
                "bp_in": data.get("bp_in"),  # number of input base-pairs
                "bp_out": data.get("bp_out"),  # number of output base-pairs
                "bp_out_proportion": data.get("bp_out_proportion"),
                "bp_removed": data.get("bp_removed"),  # number of removed base-pairs
                "bp_removed_proportion": data.get("bp_removed_proportion"),
                # "time" : data.get("time"), #runtime
                # "seqs_per_second" : data.get("seqs_per_second"),
                # "bp_per_second" : data.get("bp_per_second")
            }

        # check, if the dict is empty (if so, it either failed to parse or didnt match any deacon reports)
        if len(self.deacon_data) == 0:
            log.warning("no deacon reports found")

        # General Stats Table
        self.deacon_general_stats_table(self.deacon_data)

        # create table in the report
        self.add_section(
            name="Deacon statistics",
            anchor="deacon_stats",
            description="Statistics parsed from JSON reports",
            plot=table.plot(
                self.deacon_data,
                {
                    # "k" : {"title" : "k-mer length"},
                    # "w" : {"title" : "window size"},
                    # "abs_threshold" : {"title" : "absolute threshold"},
                    # "rel_threshold" : {"title" : "relative threshold"},
                    # "prefix_length" : {"title" : "prefix length"},
                    "deplete": {"title": "deplete"},
                    # "rename" : {"title" : "renamed seqs"},
                    "seqs_in": {"title": "Sequences in"},
                    "seqs_out": {"title": "Sequences out"},
                    "seqs_removed": {"title": "Sequences removed"},
                    "seqs_removed_proportion": {
                        "title": "Sequences removed (%)",
                        "format": "{:.2%}",
                    },
                    "bp_in": {"title": "bp in"},
                    "bp_out": {"title": "bp out"},
                    "bp_removed": {"title": "bp removed"},
                    "bp_removed_proportion": {
                        "title": "bp removed (%)",
                        "format": "{:.2%}",
                    },
                    # "time": {"title" : "time (in s)"},
                    # "seqs_per_second" : {"title" : "reads/s"},
                    # "bp_per_second" : {"title" : "bp/s"}
                },
                {
                    "id": "Deacon_statistics",
                    "title": "Deacon report statistics",
                },
            ),
        )

        # add plot in report for: Sequences removed and kept in absolute and percentage number

        # initialize an empty dict to hold data for the barplot
        plot_data = {}
        # initialize an empty list to contain the category names for the bargraph
        cats: List[str] = list()

        # iterate through all samples and their statistics
        for sample, stats in self.deacon_data.items():
            removed = stats.get("seqs_removed")  # number of removed sequences
            seqs_out = stats.get("seqs_out")  # number of kept sequences

            # only include the sample, if both metrics exist
            if removed is not None and seqs_out is not None:
                plot_data[sample] = {
                    "removed": float(removed),
                    "kept": float(seqs_out),
                    # convert both to float for plotting consistency
                }

        if plot_data:  # create if at least one sample has data
            # add category names
            cats.append("removed")
            cats.append("kept")

            # configuration of plot settings
            pconfig_plot = {
                "id": "Deacon_Sequences_removed_vs_kept_plot",
                "title": "Removed and kept sequences",
                "ylab": "Number of Sequences",
                "cpswitch": True,  # switch between absolute and percentage number
                "stacking": "relative",  # stacked bars
            }

            # add a new section in MultiQC report
            self.add_section(
                name="Removed and kept sequences",
                anchor="Deacon_Sequences_removed_vs_kept_section",
                description="Number of removed and kept sequences",
                plot=bargraph.plot(plot_data, cats, pconfig_plot),  # generates a barplot
            )

    def deacon_general_stats_table(self, data_by_sample):
        """Take the parsed stats from the fastp report and add it to the
        General Statistics table at the top of the report"""

        self.general_stats_addcols(
            data_by_sample,
            headers={
                "seqs_removed": {
                    "title": "Sequences removed",
                    "description": "Sequences removed while filtering",
                },
                "seqs_out": {
                    "title": "Sequences kept",
                    "description": "Number of sequences kept after filtering",
                },
            },
        )
