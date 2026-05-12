import json
import logging
import re

from typing import List

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics from the JSON log file (--summary or -s options). It parses and reports
    the number of input, output and removed sequences, and displays them in the General Stats table as well as a bar plot.
    """

    def __init__(self):
        super().__init__(
            name="Deacon",
            anchor="deacon",
            info="Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching.",
            href="https://github.com/bede/deacon",
            doi="https://doi.org/10.1101/2025.06.09.658732",
        )

        self.deacon_data = dict()

        for f in self.find_log_files("deacon", filehandles=True):
            try:
                data = json.load(f["f"])
            except json.JSONDecodeError as e:
                log.warning(f"Failed to parse JSON from {f['fn']}: {e}")
                continue

            version = data.get("version", "")
            if not version.startswith("deacon"):
                log.debug(f"Skipping {f['fn']}: not a Deacon report")
                continue

            self.add_data_source(f)

            try:
                self.deacon_data[f["s_name"]] = {
                    # "k" : data.get("k"), # k-mere length
                    # "w" : data.get("w"), # window size
                    # "abs_threshold" : data.get("abs_threshold"), # number of k-mers that count as hit
                    # "rel_threshold" : data.get("rel_threshold"),
                    # "prefix_length" : data.get("prefix_length"), # prefix length used for hashing/filtering
                    "deplete": "True" if data["deplete"] else "False",
                    # "rename" : data.get("rename"), # boolean, renamed seqs after export (True/False)
                    "seqs_in": data["seqs_in"],
                    "seqs_out": data["seqs_out"],
                    "seqs_out_proportion": data["seqs_out_proportion"],
                    "seqs_removed": data["seqs_removed"],
                    "seqs_removed_proportion": data["seqs_removed_proportion"],
                    "bp_in": data["bp_in"],
                    "bp_out": data["bp_out"],
                    "bp_out_proportion": data["bp_out_proportion"],
                    "bp_removed": data["bp_removed"],
                    "bp_removed_proportion": data["bp_removed_proportion"],
                    # "time" : data.get("time"), #runtime
                    # "seqs_per_second" : data.get("seqs_per_second"),
                    # "bp_per_second" : data.get("bp_per_second")
                }
            except KeyError as e:
                raise KeyError(f"Missing required field {e} in Deacon report {f['fn']}") from e

            self.add_software_version(re.sub(r"^\D*", "", version), sample=f["s_name"])

        self.deacon_data = self.ignore_samples(self.deacon_data)

        if len(self.deacon_data) == 0:
            raise ModuleNoSamplesFound

        # General Stats Table
        self.general_stats_addcols(
            self.deacon_data,
            headers={
                "seqs_removed": {
                    "title": "Sequences removed",
                    "description": "Sequences removed while filtering",
                    "scale": "Reds",
                    "shared_key": "read_count",
                },
                "seqs_out": {
                    "title": "Sequences kept",
                    "description": "Number of sequences kept after filtering",
                    "scale": "Greens",
                    "shared_key": "read_count",
                },
            },
        )

        # Detailed stats table
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
                    "deplete": {"title": "Deplete"},
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
                    "id": "deacon_statistics_table",
                    "title": "Deacon report statistics",
                },
            ),
        )

        # Bar plot: removed vs kept sequences
        plot_data = {}
        cats: List[str] = []

        for sample, stats in self.deacon_data.items():
            removed = stats.get("seqs_removed")
            seqs_out = stats.get("seqs_out")
            if removed is not None and seqs_out is not None:
                plot_data[sample] = {
                    "removed": float(removed),
                    "kept": float(seqs_out),
                }

        if plot_data:
            cats.append("removed")
            cats.append("kept")

            pconfig_plot = {
                "id": "deacon_sequences_removed_vs_kept_plot",
                "title": "Removed and kept sequences",
                "ylab": "Number of Sequences",
                "cpswitch": True,
                "stacking": "relative",
            }

            self.add_section(
                name="Removed and kept sequences",
                anchor="deacon_sequences_removed_vs_kept_section",
                description="Number of removed and kept sequences",
                plot=bargraph.plot(plot_data, cats, pconfig_plot),
            )

        # Write parsed data to file (must be after all sections)
        self.write_data_file(self.deacon_data, "multiqc_deacon")
