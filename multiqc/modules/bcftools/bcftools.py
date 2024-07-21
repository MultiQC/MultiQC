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
    [your configuration](http://multiqc.info/docs/#configuring-multiqc):

    ```yaml
    bcftools:
      collapse_complementary_changes: true
    ```

    MultiQC will sum up all complementary changes and show only `A>*` and `C>*` substitutions
    in the resulting plot.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bcftools",
            anchor="bcftools",
            target="Bcftools",
            href="https://samtools.github.io/bcftools/",
            info="Utilities for variant calling and manipulating VCFs and BCFs.",
            doi="10.1093/gigascience/giab008",
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
