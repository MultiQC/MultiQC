""" MultiQC module to parse output from Sentieon-dnaseq """


import logging
import os
import re
from typing import Dict, Optional

from multiqc.modules.picard import MultiqcModule as PicardModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(PicardModule):
    """
    Module for collecting QC stats from Sentieon-dnaseq, all of them are produced
    by Picard tools internally, though containing different headers compared to
    the original Picard tools output.
    """

    def __init__(self):
        # Inherit from Picard MultiqcModule rather than base MultiqcModule
        super().__init__(
            name="Sentieon",
            anchor="sentieon",
            href="https://www.sentieon.com/products/",
            info="contains a suite of QC tools. The ones represented in this module are analogues of Picard QC metrics.",
            # Can't find a DOI // doi=
        )

    def is_line_right_before_table(self, line: str) -> bool:
        """
        Picard logs from different samples can be concatenated together, so the module
        needs to know a marker to find where new sample information starts.

        Sentieon uses Picard tools, but adds its own header.
        """
        return line.startswith("#SentieonCommandLine:")

    def extract_sample_name(self, line: str, f: Dict) -> Optional[str]:
        """
        A file can be concatenated from multiple samples, so we can't just extract
        the sample name from the file name, and need a way to find sample name in the
        header. The copy of the originally used command is the best bet, as it's
        usually logged by Picard.

        Sentieon uses Picard tools, but adds its own header.
        """
        if line.startswith("#SentieonCommandLine:") and " --algo " in line and " -i " in line:
            # Pull sample name from the input file name, recorded in the command line:
            fn_search = re.search(r" -i\s+(\[?\S+\]?)", line, flags=re.IGNORECASE)
            if fn_search:
                s_name = os.path.basename(fn_search.group(1).strip("[]"))
                s_name = self.clean_s_name(s_name, f)
                return s_name
        return None
