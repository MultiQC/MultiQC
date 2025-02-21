import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .err_spl import parse_glimpse_err_spl
from .err_grp import parse_glimpse_err_grp

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The supported files are generated from the `GLIMPSE2_concordance` command. The following files are supported:

    - `*.error.spl.txt`
    - `*.error.grp.txt`
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GLIMPSE",
            anchor="glimpse",
            target="GLIMPSE",
            href="https://odelaneau.github.io/GLIMPSE/",
            info="Low-coverage whole genome sequencing imputation",
            extra="""
            The program `GLIMPSE2` is based on the GLIMPSE model and designed for reference panels containing
            hundreds of thousands of reference samples, with a special focus on rare variants.
        
            The concordance rates values are displayed in a scatter plot, with the option to switch between
            the different concordance metrics.
            """,
            doi="10.1101/2022.11.28.518213 ",
        )

        # Call submodule functions
        n_reports_found = 0
        n_reports_found += parse_glimpse_err_spl(self)
        n_reports_found += parse_glimpse_err_grp(self)

        # Exit if we didn't find anything
        if n_reports_found == 0:
            raise ModuleNoSamplesFound
