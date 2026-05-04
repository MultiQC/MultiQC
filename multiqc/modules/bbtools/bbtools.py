"""MultiQC module to parse output from BBTools suite"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .aqhist import parse_bbtools_aqhist
from .bbsplit import parse_bbtools_bbsplit
from .bhist import parse_bbtools_bhist
from .bqhist import parse_bbtools_bqhist
from .covhist import parse_bbtools_covhist
from .duk import parse_bbtools_duk
from .ehist import parse_bbtools_ehist
from .gchist import parse_bbtools_gchist
from .idhist import parse_bbtools_idhist
from .ihist import parse_bbtools_ihist
from .indelhist import parse_bbtools_indelhist
from .lhist import parse_bbtools_lhist
from .mhist import parse_bbtools_mhist
from .qahist import parse_bbtools_qahist
from .qchist import parse_bbtools_qchist
from .qhist import parse_bbtools_qhist
from .stats import parse_bbtools_stats

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA
    and RNA sequence data. The module parses outputs from various BBTools programs.

    Note: This module was previously split into separate `bbmap` and `bbduk` modules. These have
    been unified into a single `bbtools` module. If you were previously using module-specific
    configuration for `bbmap` or `bbduk`, you may need to update your configuration to use `bbtools`.

    Supported commands:

    - `bbduk` - Filtering and trimming statistics from stdout logs
    - `bbsplit` - Statistics on how many reads mapped to which reference genome
    - `stats` - BBDuk filtering statistics (adapter/contaminant matching)
    - `aqhist` - Histogram of average read quality
    - `bhist` - Base composition histogram by position
    - `bqhist` - Quality histogram designed for box plots
    - `covhist` - Histogram of coverage depth levels
    - `ehist` - Errors-per-read histogram
    - `gchist` - Read GC content histogram
    - `idhist` - Histogram of read identity percentages
    - `ihist` - Insert size histogram
    - `indelhist` - Indel length histogram
    - `lhist` - Read length histogram
    - `mhist` - Match, substitution, deletion, and insertion rates by position
    - `qahist` - Quality accuracy histogram
    - `qchist` - Count of bases with each quality value
    - `qhist` - Quality histogram by position
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="BBTools",
            anchor="bbtools",
            href="http://jgi.doe.gov/data-and-tools/bbtools/",
            info="Fast, multithreaded bioinformatics tools for DNA/RNA sequence analysis.",
            # One publication, but only for the merge tool:
            # doi="10.1371/journal.pone.0185056",
        )

        n = dict()

        # Call submodule functions (ordered alphabetically by section name)
        # BBDuk sections
        n["duk"] = parse_bbtools_duk(self)
        if n["duk"] > 0:
            log.info(f"Found {n['duk']} BBDuk reports")

        n["stats"] = parse_bbtools_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        # BBSplit sections
        n["bbsplit"] = parse_bbtools_bbsplit(self)
        if n["bbsplit"] > 0:
            log.info(f"Found {n['bbsplit']} BBSplit reports")

        # BBMap sections
        n["bhist"] = parse_bbtools_bhist(self)
        if n["bhist"] > 0:
            log.info(f"Found {n['bhist']} bhist reports")

        n["bqhist"] = parse_bbtools_bqhist(self)
        if n["bqhist"] > 0:
            log.info(f"Found {n['bqhist']} bqhist reports")

        n["qchist"] = parse_bbtools_qchist(self)
        if n["qchist"] > 0:
            log.info(f"Found {n['qchist']} qchist reports")

        n["covhist"] = parse_bbtools_covhist(self)
        if n["covhist"] > 0:
            log.info(f"Found {n['covhist']} covhist reports")

        n["ehist"] = parse_bbtools_ehist(self)
        if n["ehist"] > 0:
            log.info(f"Found {n['ehist']} ehist reports")

        n["gchist"] = parse_bbtools_gchist(self)
        if n["gchist"] > 0:
            log.info(f"Found {n['gchist']} gchist reports")

        n["idhist"] = parse_bbtools_idhist(self)
        if n["idhist"] > 0:
            log.info(f"Found {n['idhist']} idhist reports")

        n["indelhist"] = parse_bbtools_indelhist(self)
        if n["indelhist"] > 0:
            log.info(f"Found {n['indelhist']} indelhist reports")

        n["ihist"] = parse_bbtools_ihist(self)
        if n["ihist"] > 0:
            log.info(f"Found {n['ihist']} ihist reports")

        n["mhist"] = parse_bbtools_mhist(self)
        if n["mhist"] > 0:
            log.info(f"Found {n['mhist']} mhist reports")

        n["qahist"] = parse_bbtools_qahist(self)
        if n["qahist"] > 0:
            log.info(f"Found {n['qahist']} qahist reports")

        n["qhist"] = parse_bbtools_qhist(self)
        if n["qhist"] > 0:
            log.info(f"Found {n['qhist']} qhist reports")

        n["lhist"] = parse_bbtools_lhist(self)
        if n["lhist"] > 0:
            log.info(f"Found {n['lhist']} lhist reports")

        n["aqhist"] = parse_bbtools_aqhist(self)
        if n["aqhist"] > 0:
            log.info(f"Found {n['aqhist']} aqhist reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
