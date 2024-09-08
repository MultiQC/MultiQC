import logging
from typing import Dict, List, Union

from multiqc.modules.kraken import MultiqcModule as KrakenModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(KrakenModule):
    """
    This module works with Bracken output files that resemble Kraken reports. They look like the following:

    ```tsv
    100.00	1188381	0	R	1	root
    100.00	1188380	0	R1	131567	  cellular organisms
    100.00	1188328	0	D	2	    Bacteria
    99.99	1188317	0	P	1224	      Proteobacteria
    99.99	1188308	0	C	1236	        Gammaproteobacteria
    99.99	1188285	0	O	91347	          Enterobacterales
    99.98	1188147	0	F	543	            Enterobacteriaceae
    ```

    The main assumption to tell Bracken reports from Kraken is that the former don't have
    an "unassigned" category at the head of the file, and instead start with "root".

    A bar graph is generated identical to that of Kraken.
    """

    def __init__(self):
        super().__init__(
            name="Bracken",
            anchor="bracken",
            href="https://ccb.jhu.edu/software/bracken/",
            info="Computes the abundance of species in DNA sequences from a metagenomics sample.",
            doi="10.7717/peerj-cs.104",
            sp_key="bracken",
        )
