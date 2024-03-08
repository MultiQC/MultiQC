""" MultiQC module to parse output from bracken """

import logging

from multiqc.modules.kraken import MultiqcModule as KrakenModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(KrakenModule):
    """Bracken module"""

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="Bracken",
            anchor="bracken",
            href="https://ccb.jhu.edu/software/bracken/",
            info="is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample.",
            doi="10.7717/peerj-cs.104",
            sp_key="bracken",
        )
