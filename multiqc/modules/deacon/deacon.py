import json
import logging
import re

from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    [Deacon](https://github.com/bede/deacon) filters DNA sequences in FASTA/Q files and streams
    using SIMD-accelerated minimizer comparison against an indexed query. It can either keep
    matching sequences (search mode, default) or remove them (depletion mode, `-d` / `--deplete`).
    Built with panhuman host depletion in mind but useful for searching large sequence
    collections.

    This module parses the JSON summary log written by `deacon filter` when called with
    `-s` / `--summary`, and reports the number of input, kept, and removed sequences and base
    pairs alongside the filter mode.

    ## Generating compatible output

    The module looks for the JSON summary file produced by the `--summary` (`-s`) option:

    ```bash
    # Search mode: keep matching reads
    deacon filter index.idx reads.fq.gz -o matches.fq.gz -s summary.json

    # Depletion mode: remove matching reads (e.g. host depletion)
    deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o depleted.fq.gz -s summary.json
    ```

    ## Interpreting results

    `seqs_removed` and `bp_removed` count sequences and base pairs that matched the indexed
    query and were therefore filtered out. In depletion mode (`Deplete = True`) these are
    typically host reads to discard; in search mode they are non-target reads. The `Deplete`
    column distinguishes the two modes so the same report can mix search and depletion
    samples.
    """

    def __init__(self):
        super().__init__(
            name="Deacon",
            anchor="deacon",
            info="Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching.",
            href="https://github.com/bede/deacon",
            doi="https://doi.org/10.1101/2025.06.09.658732",
        )

        self.deacon_data: Dict[str, Dict] = {}

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
                    "title": f"Seqs removed ({config.read_count_prefix})",
                    "description": f"Number of sequences removed by the filter ({config.read_count_desc})",
                    "scale": "Reds",
                    "shared_key": "read_count",
                    "modify": lambda x: float(x) * config.read_count_multiplier,
                    "min": 0,
                    "tt_decimals": 2,
                },
                "seqs_removed_proportion": {
                    "title": "% Removed",
                    "description": "Percentage of sequences removed by the filter",
                    "scale": "OrRd",
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "modify": lambda x: float(x) * 100.0,
                    "tt_decimals": 2,
                    "hidden": True,
                },
                "seqs_out": {
                    "title": f"Seqs kept ({config.read_count_prefix})",
                    "description": f"Number of sequences passing the filter ({config.read_count_desc})",
                    "scale": "Greens",
                    "shared_key": "read_count",
                    "modify": lambda x: float(x) * config.read_count_multiplier,
                    "min": 0,
                    "tt_decimals": 2,
                },
            },
            namespace="Deacon",
        )

        # Detailed stats table
        self.add_section(
            name="Deacon statistics",
            anchor="deacon_stats",
            description="Sequence and base counts before and after filtering, parsed from Deacon JSON summaries.",
            helptext="""
            **Columns:**

            * **Deplete** — `True` if the sample was filtered with `--deplete` (`-d`), in which
              case matching sequences were removed (e.g. host depletion). `False` means search
              mode, where matching sequences were kept.
            * **Sequences in / out / removed** — Counts of input sequences, sequences passing
              the filter, and sequences filtered out, respectively. For paired-end data each
              read pair is counted as a single sequence.
            * **Sequences removed (%)** — Proportion of input sequences that were filtered out.
            * **Bases in / out / removed** — Equivalent counts in base pairs.
            * **Bases removed (%)** — Proportion of input base pairs that were filtered out.

            In search mode, `Sequences removed` represents reads that did not match the indexed
            query. In depletion mode, `Sequences removed` represents reads that did match (and
            were therefore discarded as e.g. host contamination).
            """,
            plot=table.plot(
                self.deacon_data,
                {
                    # "k" : {"title" : "k-mer length"},
                    # "w" : {"title" : "window size"},
                    # "abs_threshold" : {"title" : "absolute threshold"},
                    # "rel_threshold" : {"title" : "relative threshold"},
                    # "prefix_length" : {"title" : "prefix length"},
                    "deplete": {
                        "title": "Deplete",
                        "description": "Filter mode: True for depletion (remove matches), False for search (keep matches)",
                        "scale": False,
                    },
                    # "rename" : {"title" : "renamed seqs"},
                    "seqs_in": {
                        "title": f"Seqs in ({config.read_count_prefix})",
                        "description": f"Total input sequences ({config.read_count_desc})",
                        "scale": "Blues",
                        "shared_key": "read_count",
                        "modify": lambda x: float(x) * config.read_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "seqs_out": {
                        "title": f"Seqs kept ({config.read_count_prefix})",
                        "description": f"Sequences passing the filter ({config.read_count_desc})",
                        "scale": "Greens",
                        "shared_key": "read_count",
                        "modify": lambda x: float(x) * config.read_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "seqs_removed": {
                        "title": f"Seqs removed ({config.read_count_prefix})",
                        "description": f"Sequences removed by the filter ({config.read_count_desc})",
                        "scale": "Reds",
                        "shared_key": "read_count",
                        "modify": lambda x: float(x) * config.read_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "seqs_removed_proportion": {
                        "title": "Seqs removed",
                        "description": "Percentage of input sequences removed by the filter",
                        "scale": "OrRd",
                        "min": 0,
                        "max": 100,
                        "suffix": "%",
                        "modify": lambda x: float(x) * 100.0,
                        "tt_decimals": 2,
                    },
                    "bp_in": {
                        "title": f"Bases in ({config.base_count_prefix})",
                        "description": f"Total input base pairs ({config.base_count_desc})",
                        "scale": "Blues",
                        "shared_key": "base_count",
                        "modify": lambda x: float(x) * config.base_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "bp_out": {
                        "title": f"Bases kept ({config.base_count_prefix})",
                        "description": f"Base pairs passing the filter ({config.base_count_desc})",
                        "scale": "Greens",
                        "shared_key": "base_count",
                        "modify": lambda x: float(x) * config.base_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "bp_removed": {
                        "title": f"Bases removed ({config.base_count_prefix})",
                        "description": f"Base pairs removed by the filter ({config.base_count_desc})",
                        "scale": "Reds",
                        "shared_key": "base_count",
                        "modify": lambda x: float(x) * config.base_count_multiplier,
                        "min": 0,
                        "tt_decimals": 2,
                    },
                    "bp_removed_proportion": {
                        "title": "Bases removed",
                        "description": "Percentage of input base pairs removed by the filter",
                        "scale": "OrRd",
                        "min": 0,
                        "max": 100,
                        "suffix": "%",
                        "modify": lambda x: float(x) * 100.0,
                        "tt_decimals": 2,
                    },
                    # "time": {"title" : "time (in s)"},
                    # "seqs_per_second" : {"title" : "reads/s"},
                    # "bp_per_second" : {"title" : "bp/s"}
                },
                {
                    "id": "deacon_statistics_table",
                    "title": "Deacon: Filter statistics",
                },
            ),
        )

        # Bar plot: removed vs kept sequences
        plot_data: Dict[str, Dict[str, float]] = {}
        for sample, stats in self.deacon_data.items():
            removed = stats.get("seqs_removed")
            seqs_out = stats.get("seqs_out")
            if removed is not None and seqs_out is not None:
                plot_data[sample] = {
                    "seqs_out": float(seqs_out),
                    "seqs_removed": float(removed),
                }

        if plot_data:
            cats = {
                "seqs_out": {"name": "Sequences kept", "color": "#7cb5ec"},
                "seqs_removed": {"name": "Sequences removed", "color": "#f15c80"},
            }

            pconfig_plot = {
                "id": "deacon_sequences_removed_vs_kept_plot",
                "title": "Deacon: Sequences kept vs removed",
                "ylab": "Number of sequences",
                "cpswitch_counts_label": "Number of sequences",
                "cpswitch_percent_label": "Percentage of sequences",
            }

            self.add_section(
                name="Sequences kept vs removed",
                anchor="deacon_sequences_removed_vs_kept_section",
                description="Sequences passing and failing the Deacon filter, per sample.",
                helptext="""
                Each bar shows the number of sequences kept and removed by Deacon for a
                given sample. Use the toggle above the plot to switch between absolute
                counts and percentages.

                Whether the kept or removed bar is the desired outcome depends on the
                filter mode: in depletion mode (e.g. host depletion) you typically want
                most reads to be **kept**, while in search mode you typically expect only
                a small **kept** fraction matching the indexed query.
                """,
                plot=bargraph.plot(plot_data, cats, pconfig_plot),
            )

        # Write parsed data to file (must be after all sections)
        self.write_data_file(self.deacon_data, "multiqc_deacon")
