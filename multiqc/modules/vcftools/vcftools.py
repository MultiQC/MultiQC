#!/usr/bin/env python

""" MultiQC module to parse output from vcftools """

from __future__ import print_function
import logging
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from .relatedness2 import Relatedness2Mixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, Relatedness2Mixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='VCFTools',
            anchor='vcftools',
            href='https://vcftools.github.io',
            info='is a program for working with and reporting on VCF files.'
        )

        n = dict()
        n['relatedness2'] = self.parse_relatedness2()

        if sum(n.values()) == 0:
            log.debug('Could not find any reports in %s', config.analysis_dir)
            raise UserWarning
