import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .compare import CompareMixin
from .gather import GatherMixin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, CompareMixin, GatherMixin):
    """
    The module can summarise data from the following sourmash output files
    (descriptions from command line help output):

    - `sourmash compare`
      - create a similarity matrix comparing many samples.
    - `sourmash gather`
      - search a metagenome signature against databases.

    Additional information on sourmash and its outputs is available on
    the [sourmash documentation website](https://sourmash.readthedocs.io/en/latest/).

    `sourmash gather` is modelled after the Kraken module, and builds a bar graph that
    shows the coverage of top-5 genomes covered most by all samples. The number of top
    genomes can be customized in the config file:

    ```yaml
    sourmash:
      gather:
        top_n: 5
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sourmash",
            anchor="sourmash",
            href="https://github.com/sourmash-bio/sourmash",
            info="Quickly searches, compares, and analyzes genomic and metagenomic data sets.",
            doi="10.21105/joss.00027",
        )

        n = dict()
        n["compare"] = self.parse_compare()
        if n["compare"] > 0:
            log.info(f"Found {n['compare']} compare results")

        n["gather"] = self.parse_gather()
        if n["gather"] > 0:
            log.info(f"Found {n['gather']} gather results")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
