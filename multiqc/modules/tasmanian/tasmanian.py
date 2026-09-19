"""MultiQC module to parse output from Tasmanian"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .inconsistencies import parse_tasmanian_inconsistencies
from .mismatch import parse_tasmanian_mismatch
from .variants import parse_tasmanian_variants

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Tasmanian-Mismatch characterizes sequencing artifacts by counting mismatches between aligned
    reads and the reference, by position in the read or in the sequenced fragment. Damage from
    FFPE, oxidation, deamination, or other damage may appear as position-dependent
    substitution patterns that are distinct from true variation.

    Supported commands:

    - `tasmanian-mismatch`
    - `tasmanian-diagnostics`

    #### tasmanian-mismatch

    The module parses the tab-delimited mismatch count table written with `--position-mode insert`
    (the default). Tables written with `--position-mode read` or `--normalize` are not supported.
    It plots the percentage of each substitution class at every normalized fragment position, as
    a line plot and as heatmaps that compare many libraries, and adds the overall mismatch rate
    to the general statistics table. Counts are summed over the `reference_order` column, and
    reads 1 and 2 are merged into a single profile.

    ```bash
    tasmanian-mismatch sample.bam reference.fa --position-mode insert -o sample.mismatch.tsv
    ```

    #### tasmanian-diagnostics

    The two tables written by `tasmanian-diagnostics` are parsed when present:

    - `potential_variants.tsv`: genomic sites with recurrent mismatches, summarized as sites per substitution class
    - `read_pair_inconsistencies.tsv`: positions where overlapping mates disagree, summarized by discordance type

    These files have the same name for every sample by default, so give them sample-specific
    names or place them in separate directories (and run MultiQC with `-d`).

    ```bash
    tasmanian-diagnostics sample.bam reference.fa \\
      --variants-output sample.potential_variants.tsv \\
      --inconsistencies-output sample.read_pair_inconsistencies.tsv
    ```
    """

    def __init__(self):
        super().__init__(
            name="Tasmanian",
            anchor="tasmanian",
            href="https://github.com/nebiolabs/tasmanian-mismatch",
            info="Characterizes sequencing artifacts from mismatch patterns by normalized fragment position.",
            # No publication for the Rust rewrite yet
            doi=None,
            license="GNU Affero General Public License v3.0",
            license_url="https://github.com/nebiolabs/tasmanian-mismatch/blob/master/LICENCE.txt",
        )

        n = {
            "mismatch": parse_tasmanian_mismatch(self),
            "variants": parse_tasmanian_variants(self),
            "inconsistencies": parse_tasmanian_inconsistencies(self),
        }
        for name, count in n.items():
            if count > 0:
                log.info(f"Found {count} {name} reports")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
