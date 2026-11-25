import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .counts import parse_dexseq_counts

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    DEXSeq tests RNA-seq data for differential exon usage. The module parses the exon bin
    counts written by `dexseq_count.py`.

    `dexseq_count.py` ships with the DEXSeq R package. It counts the reads overlapping each exon
    bin (exonic part) of an annotation flattened by `dexseq_prepare_annotation.py`, and writes one
    file per sample, ending with counters of the reads it could not assign:

    ```bash
    python dexseq_count.py --paired yes --format bam annotation.gff sample.bam sample.txt
    ```

    Sample names are taken from the file names. The module shows the number of exon bins detected
    in each sample, and adds the unassigned reads to the General Statistics table (hidden by
    default). A read overlapping several exon bins of one gene is counted once for each bin, and
    the script does not report how many reads it counted, so the module does not give a
    percentage of assigned reads.
    """

    def __init__(self):
        super().__init__(
            name="DEXSeq",
            anchor="dexseq",
            href="https://bioconductor.org/packages/DEXSeq/",
            info="Tests RNA-seq data for differential exon usage.",
            doi="10.1101/gr.133744.111",
            license="GNU General Public License v3.0 or later",
            license_url="https://www.gnu.org/licenses/gpl-3.0.html",
        )

        n_counts = parse_dexseq_counts(self)
        if n_counts == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {n_counts} exon bin count files")
