#!/usr/bin/env python

""" MultiQC module to parse output from Lima """

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
          name='Freyja',
          anchor='freyja',
          href="https://github.com/andersen-lab/Freyja",
          info="Recover relative lineage abundances from mixed SARS-CoV-2 samples."
        )