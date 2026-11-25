import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .read_outcomes import parse_rmats_read_outcomes
from .summary import parse_rmats_summary

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    rMATS-turbo detects differential alternative splicing events between two groups of
    RNA-seq samples. It reports five event types: skipped exon (SE), alternative 5' and
    3' splice sites (A5SS, A3SS), mutually exclusive exons (MXE) and retained intron (RI).

    The module parses two rMATS-turbo outputs (v4.1.1 and later).

    #### Read outcomes

    The `prep` step writes `[datetime]_read_outcomes_by_bam.txt` to the `--tmp`
    directory. For each input BAM it records how many reads were used and why the
    others were discarded. The sample name is taken from the BAM file name in the file,
    so a file covering several BAMs gives one sample per BAM. When the BAM is named
    `Aligned.sortedByCoord.out.bam` (as when rMATS maps FASTQ input itself), the name of
    the BAM's directory is used instead.

    The file name must end in `read_outcomes_by_bam.txt`. rMATS-turbo's own Nextflow,
    WDL and CWL workflows rename it to `prep_<id>_read_outcomes_by_bam.txt`, which also
    matches.

    #### Summary

    `summary.txt` in the output directory counts the events tested and the events
    called significant, per event type, for both the junction count (JC) and the
    junction count plus exon count (JCEC) analyses. By default an event is significant
    when its FDR is at most 0.05. `rMATS_P/summary.py` can rewrite the summary with
    other cutoffs; with `--summary-prefix` the cutoffs are added to the file name
    (`<prefix>_FDR_0.05_IncLevelDifference_0.1.txt`) and the module shows them.

    Summary files are found by their header, so they may have any name. Each file is
    one comparison (one `--b1` versus `--b2` run), named after the file without the
    `.rmats_summary` or `_summary` suffix, for example `PC3E-GS689.rmats_summary.txt`
    becomes `PC3E-GS689`. A file named only `summary.txt` is named after its directory.
    A cutoff suffix stays in the name, so that summaries of one comparison at several
    cutoffs stay apart.

    Summary values describe comparisons, not samples, so they are not added to the
    General Statistics table.
    """

    def __init__(self):
        super().__init__(
            name="rMATS",
            anchor="rmats",
            href="https://github.com/Xinglab/rmats-turbo",
            info="Detects differential alternative splicing events from RNA-seq data.",
            doi=["10.1073/pnas.1419161111", "10.1038/s41596-023-00944-2"],
            # Free for academic and non-profit use under the 2-Clause BSD License;
            # commercial use needs a separate licence from the authors.
            license="Custom (2-Clause BSD License for academic and non-profit use)",
            license_url="https://github.com/Xinglab/rmats-turbo/blob/master/LICENSE",
        )

        n_read_outcomes = parse_rmats_read_outcomes(self)
        if n_read_outcomes > 0:
            log.info(f"Found {n_read_outcomes} BAM read outcome reports")

        n_summary = parse_rmats_summary(self)
        if n_summary > 0:
            log.info(f"Found {n_summary} summary reports")

        if n_read_outcomes + n_summary == 0:
            raise ModuleNoSamplesFound
