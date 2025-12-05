"""MultiQC module to parse output from SeqKit"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .stats import parse_seqkit_stats

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    SeqKit is a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.

    Supported commands:

    - `stats`

    The module parses output from `seqkit stats` which provides simple statistics of
    FASTA/Q files including sequence counts, total length, N50, GC content, and quality
    metrics for FASTQ files.

    #### stats

    The `seqkit stats` command produces tabular output with columns for file, format,
    type, num_seqs, sum_len, min_len, avg_len, max_len, and optionally Q1, Q2, Q3,
    sum_gap, N50, Q20(%), Q30(%), AvgQual, and GC(%) when run with the `--all` flag.

    To generate output suitable for MultiQC, run seqkit stats with the `--tabular` flag:

    ```bash
    seqkit stats --all --tabular *.fastq.gz > seqkit_stats.tsv
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="SeqKit",
            anchor="seqkit",
            href="https://bioinf.shenwei.me/seqkit/",
            info="Cross-platform and ultrafast toolkit for FASTA/Q file manipulation.",
            doi="10.1371/journal.pone.0163962",
        )

        n = dict()

        # Call submodule functions
        n["stats"] = parse_seqkit_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
