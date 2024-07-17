import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.picard import MarkDuplicates

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Currently, the biobambam2 module only processes output from the `bamsormadup` command.
    Not only that, but it cheats by using the module code from Picard/MarkDuplicates.
    The output is so similar that the code simply sets up a module with unique name and
    filename search pattern and then uses the parsing code from the Picard module.

    Apart from behind the scenes coding, this module should work in exactly the same way
    as all other MultiQC modules.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="biobambam2",
            anchor="biobambam2",
            href="https://gitlab.com/german.tischler/biobambam2",
            info="Tools for early stage alignment file processing",
            doi="10.1186/1751-0473-9-13",
        )

        n = dict()
        n["bamsormadup"] = MarkDuplicates.parse_reports(self, "biobambam2/bamsormadup")
        if len(n["bamsormadup"]) > 0:
            log.info(f"Found {len(n['bamsormadup'])} bamsormadup reports")
        else:
            raise ModuleNoSamplesFound

    # Helper functions
    @staticmethod
    def multiply_hundred(val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
