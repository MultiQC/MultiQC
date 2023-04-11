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

        # To store the summary data
        self.freyja_data = dict()

        # Parse the output files
        self.parse_summ_files()

        # Remove filtered samples
        self.freyja_data = self.ignore_samples(self.freyja_data)

        # Let MultiQC know this module found no data
        if len(self.freyja_data) == 0:
            raise UserWarning

        log.info(f"Found {len(self.freyja_data)} reports")
        self.write_data_file(self.freyja_data, "multiqc_freyja")

    def parse_summ_files(self, f):
        """
        Parse the summary file. 
        Freyja has multiple summary files, but we only need to parse the one from the demix command.
        More specifically, we only need the line that starts with "summarized".
        ...
        summarized	[('BQ.1*', 0.983), ('Omicron', 0.011), ('key', value)]
        ...
        """    
        for f in self.find_log_files('freyja'):
            s_name = self.clean_s_name(f["root"], f)

            # Read the statistics from file
            d = {}
            for line in f["f"]:
                try:
                    if line.startswith('summarized'):
                        summarized_line = line
                        summarized_line = summarized_line.strip().split('\t')[1]
                        d = eval(summarized_line) # Make sure no input file does not contain any malicious code
                        d = dict(d)
                except ValueError:
                    pass
            
            # Percentages don't always add up to 1, show a warning if this is the case
            if sum(d.values()) != 1 :
                log.warning(f"Freyja {s_name}: percentages don't sum to 1")
            
            # There is no sample name in the log, so we use the root of the
            # file as sample name (since the filename is always stats.dat
            if s_name in self.freyja_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.freyja_data[s_name] = data
            self.add_data_source(f, s_name)

    def general_stats_cols(self):
        """Add a single column displaying the most abundant variant to the General Statistics table"""
        top_variant_dict= {}
        for s_name, sub_dict in self.freyja_data.items():
            top_variant = max(sub_dict, key=sub_dict.get)
            top_variant_value = sub_dict[top_variant]
            top_variant_dict[s_name] = {
                "Top variant": top_variant,
                "Top variant %": top_variant_value
            }
        
        headers = OrderedDict()
        headers['Top variant'] = {
            'title': 'Top variant',
            'description': 'The most abundant variant in the sample',
            'scale': 'RdYlGn-rev' # Not sure if this is the best scale
        }
        headers['Top variant %'] = {
            'title': 'Top variant %',
            'description': 'The percentage of the most abundant variant in the sample',
            'max': 100,
            'min': 0,
            'scale': 'Blues',
            'suffix': '%'
        }

        self.general_stats_addcols(top_variant_dict, headers)