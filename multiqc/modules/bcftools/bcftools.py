import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from multiqc.modules.bcftools.stats import parse_bcftools_stats

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands: `stats`

    #### Collapse complementary substitutions

    In non-strand-specific data, reporting the total numbers of occurences for both changes
    in a comlementary pair - like `A>C` and `T>G` - might not bring any additional information.
    To collapse such statistics in the substitutions plot, you can add the following section into
    [your configuration](https://docs.seqera.io/multiqc/getting_started/config):

    ```yaml
    bcftools:
      collapse_complementary_changes: true
    ```

    MultiQC will sum up all complementary changes and show only `A>*` and `C>*` substitutions
    in the resulting plot.

    #### Aggregating Hom/Het counts in General Stats

    When `bcftools stats` is run with `-s` or `-S`, the PSC section contains homozygous
    and heterozygous counts for each sample in the input file. MultiQC currently represents
    each bcftools set as one row in General Stats, so these per-sample values are aggregated
    for the `Hom` and `Het` columns. The default is the rounded mean; this can be changed to
    the rounded median or sum:

    ```yaml
    bcftools:
      hom_het_aggregation: mean # or: median, sum
    ```

    `Hom` counts non-reference homozygous genotypes only, as hom-ref genotypes match the
    reference and are not variants.
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="Bcftools",
            anchor="bcftools",
            target="Bcftools",
            href="https://samtools.github.io/bcftools/",
            info="Utilities for variant calling and manipulating VCFs and BCFs.",
            doi="10.1093/gigascience/giab008",
            license="MIT License",
            license_url="https://github.com/samtools/bcftools/blob/master/LICENSE",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers: Dict = dict()
        self.general_stats_data: Dict = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = parse_bcftools_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
