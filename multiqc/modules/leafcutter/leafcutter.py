from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .clustering import parse_leafcutter_clustering
from .differential_splicing import parse_leafcutter_cluster_significance


class MultiqcModule(BaseMultiqcModule):
    """
    LeafCutter quantifies RNA splicing from the introns spanned by split reads, grouping
    introns that share splice sites into clusters, and tests clusters for differential
    splicing between two groups of samples.

    The module supports output from both the original R implementation
    ([davidaknowles/leafcutter](https://github.com/davidaknowles/leafcutter)) and the Python
    reimplementation of the differential splicing step
    ([leafcutter2/leafcutter-ds](https://github.com/leafcutter2/leafcutter-ds), packaged on
    Bioconda as `leafcutter` 2.x). The following files are parsed:

    - `*_perind_numers.counts.gz`, from `leafcutter_cluster_regtools.py` or `leafcutter_cluster.py`:
      per-sample junction read counts for each clustered intron. One file holds every sample,
      named in the header after its junction file (without the `.junc` extension).
      The module reports the introns with at least one read, and the junction reads in clusters,
      for each sample. These files are always gzip-compressed and are matched by filename only.
    - `*_cluster_significance.txt`, from `leafcutter_ds.R` or `leafcutter-ds`: the test status
      and adjusted p-value of each cluster. Each file is one comparison, named after the
      `--output_prefix` given to LeafCutter.

    The comparisons are not samples, so the differential splicing results are not added to the
    General Statistics table.

    #### Configuration

    Set this option under `leafcutter` in a MultiQC config file:

    - `padj_cutoff` (default `0.05`): a tested cluster is significant when its adjusted
      p-value is below this value. The default is the FDR threshold of LeafCutter's
      `prepare_results.R`.

    ```yaml
    leafcutter:
      padj_cutoff: 0.05
    ```

    Example commands:

    ```bash
    python leafcutter_cluster_regtools.py -j juncfiles.txt -o mystudy
    leafcutter_ds.R --num_threads 4 --exon_file exons.txt.gz -o treated_vs_control \\
        mystudy_perind_numers.counts.gz groups.txt
    ```
    """

    def __init__(self):
        super().__init__(
            name="LeafCutter",
            anchor="leafcutter",
            href=["https://davidaknowles.github.io/leafcutter/", "https://github.com/leafcutter2/leafcutter-ds"],
            info="Annotation-free quantification and differential analysis of RNA splicing from intron clusters.",
            doi="10.1038/s41588-017-0004-9",
            license="Apache License 2.0",
            license_url="https://github.com/davidaknowles/leafcutter/blob/master/LICENSE",
        )

        n_reports_found = 0
        n_reports_found += parse_leafcutter_clustering(self)
        n_reports_found += parse_leafcutter_cluster_significance(self)

        if n_reports_found == 0:
            raise ModuleNoSamplesFound
