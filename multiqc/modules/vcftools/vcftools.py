#!/usr/bin/env python

""" MultiQC module to parse output from vcftools """

from __future__ import print_function
import logging
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from .relatedness2 import Relatedness2Mixin
from .tstv_by_count import TsTvByCountMixin
from .tstv_by_qual import TsTvByQualMixin
from .tstv_summary import TsTvSummaryMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, Relatedness2Mixin, TsTvByCountMixin, TsTvByQualMixin, TsTvSummaryMixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='VCFTools',
            anchor='vcftools',
            href='https://vcftools.github.io',
            info='is a program for working with and reporting on VCF files.'
        )

        n = dict()
        n['relatedness2'] = self.parse_relatedness2()

        n['tstv_by_count'] = self.parse_tstv_by_count()
        if n['tstv_by_count'] > 0:
            log.info("Found {} TsTv.count reports".format(n['tstv_by_count']))

        n['tstv_by_qual'] = self.parse_tstv_by_qual()
        if n['tstv_by_qual'] > 0:
            log.info("Found {} TsTv.qual reports".format(n['tstv_by_qual']))

        n['tstv_summary'] = self.parse_tstv_summary()
        if n['tstv_summary'] > 0:
            log.info("Found {} TsTv.summary reports".format(n['tstv_summary']))

        if sum(n.values()) == 0:
            raise UserWarning
