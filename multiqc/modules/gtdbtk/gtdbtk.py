import logging

from multiqc.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    The module parses summary.tsv outputs from 
    """
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GTDB-Tk",
            anchor="gtdbtk",
            href="https://ecogenomics.github.io/GTDBTk/index.html",
            info="Toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.",
            doi=["10.1093/bioinformatics/btac672"],
        )
    log.info('Hello World!')