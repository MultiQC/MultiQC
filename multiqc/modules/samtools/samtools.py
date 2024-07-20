import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .stats import parse_samtools_stats
from .flagstat import parse_samtools_flagstat
from .idxstats import parse_samtools_idxstats
from .rmdup import parse_samtools_rmdup
from .coverage import parse_samtools_coverage
from .markdup import parse_samtools_markdup

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands:

    - `stats`
    - `flagstats`
    - `idxstats`
    - `rmdup`
    - `coverage`
    - `markdup`

    #### idxstats

    The `samtools idxstats` prints its results to standard out (no consistent file name) and has no header lines
    (no way to recognise from content of file). As such, `idxstats` result files must have the string `idxstat`
    somewhere in the filename.

    There are a few MultiQC config options that you can add to customise how the idxstats module works. A typical
    configuration could look as follows:

    ```yaml
    # Always include these chromosomes in the plot
    samtools_idxstats_always:
      - X
      - Y

    # Never include these chromosomes in the plot
    samtools_idxstats_ignore:
      - MT

    # Threshold where chromosomes are ignored in the plot.
    # Should be a fraction, default is 0.001 (0.1% of total)
    samtools_idxstats_fraction_cutoff: 0.001

    # Name of the X and Y chromosomes.
    # If not specified, MultiQC will search for any chromosome
    # names that look like x, y, chrx or chry (case-insensitive search)
    samtools_idxstats_xchr: myXchr
    samtools_idxstats_ychr: myYchr
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Samtools",
            anchor="samtools",
            target="samtools",
            href="http://www.htslib.org",
            info="Toolkit for interacting with BAM/CRAM files.",
            doi="10.1093/bioinformatics/btp352",
        )

        n = dict()

        # Call submodule functions
        n["stats"] = parse_samtools_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        n["flagstat"] = parse_samtools_flagstat(self)
        if n["flagstat"] > 0:
            log.info(f"Found {n['flagstat']} flagstat reports")

        n["idxstats"] = parse_samtools_idxstats(self)
        if n["idxstats"] > 0:
            log.info(f"Found {n['idxstats']} idxstats reports")

        n["rmdup"] = parse_samtools_rmdup(self)
        if n["rmdup"] > 0:
            log.info(f"Found {n['rmdup']} rmdup reports")

        n["coverage"] = parse_samtools_coverage(self)
        if n["coverage"] > 0:
            log.info(f"Found {n['coverage']} coverage reports")

        n["markdup"] = parse_samtools_markdup(self)
        if n["markdup"] > 0:
            log.info(f"Found {n['markdup']} markdup reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
