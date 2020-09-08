#!/usr/bin/env python

""" MultiQC module to parse output from CCS """

from multiqc.modules.base_module import BaseMultiqcModule


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(name='ccs', anchor='ccs',
                href='https://github.com/PacificBiosciences/ccs',
                info='Generate highly accurate single-molecule consensus '
                     'reads')
