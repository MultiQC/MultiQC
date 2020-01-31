#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module parses various
    stats produced by pairtools."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pairtools', anchor='pairtools',
        href="https://github.com/mirnylab/pairtools",
        info="pairtools is a command-line framework to process sequencing data from a Hi-C experiment.")


        # Find and load any pairtools stats summaries
        self.pairtools_stats = dict()
        for f in self.find_log_files('pairtools', filehandles=True):
            s_name = f['s_name']
            self.pairtools_stats[s_name] = self.parse_pairtools_stats(f)


        # Filter to strip out ignored sample names
        self.pairtools_stats = self.ignore_samples(self.pairtools_stats)

        if len(self.pairtools_stats) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.pairtools_stats)))
        # self.write_data_file(self.pairtools_stats, 'multiqc_pairtools')

        self.pairtools_general_stats()



    def parse_pairtools_stats(self, f):
        """ Parse a pairtools summary stats """
        # s_name = f['s_name']
        # f_name = f['fn']
        # log.info("parsing {} {} ...".format(s_name,f_name))
        #
        # just testing displying random fields from stats ...
        f_handle = f['f']
        _data = dict()
        _i = 0
        num_fields = 8
        for line in f_handle:
            _1,_2 = line.rstrip().split('\t')
            _data[_1] = int(_2)
            _i += 1
            if _i >=num_fields:
                break
        return _data


    def pairtools_general_stats(self):
        """ Add columns to General Statistics table """
        # headers = OrderedDict()
        # headers['total'] = {
        #     'title': 'total',
        #     'description': 'total number of pairs per sample',
        #     'min': 0,
        # }
        # self.general_stats_addcols(self.pairtools_stats, headers, 'pairtools')
        self.general_stats_addcols(self.pairtools_stats)
