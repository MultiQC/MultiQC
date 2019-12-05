#!/usr/bin/env python

""" MultiQC module to parse output from pycoQC """

from multiqc.modules.base_module import BaseMultiqcModule
import logging
log = logging.getLogger(__name__)




class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pycoQC', anchor='pycoqc',
        href="https://a-slide.github.io/pycoQC/",
        info="Computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data")
        log.info('Hello World!')

