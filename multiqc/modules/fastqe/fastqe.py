# https://toolshed.g2.bx.psu.edu/repository?repository_id=13576f42f394cfb6&changeset_revision=c1a33f603da2
# https://github.com/fastqe/fastqe/issues/11
import logging
import os
import re
from typing import Dict, Optional, List

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.plotly.line import LinePlotConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    **Note** - MultiQC parses the standard out from FastQE. As such, you must capture the
    FastQE stdout to a file when running to use the MultiQC module.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="FastQE",
            anchor="fastqe",
            href="https://github.com/fastqe/fastqe",
            info="FASTQ sequence quality visualisation with Emoji",
        )