""" MultiQC module to parse output from HOPS postprocessing script """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc.plots import heatmap
from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HOPS', anchor='hops',
        href="https://www.https://github.com/rhuebler/HOPS/",
        info="is an ancient DNA characteristics screening tool of output from the metagenomic aligner MALT.")

        # Find and load any HOPS post-processing JSONs
        self.hops_data = dict()

        for f in self.find_log_files('hops',filehandles=True):
            try:
                self.parseJSON(f)
            except KeyError:
                logging.warning("Error loading file {}".format(f['fn']))

        self.hops_data = self.ignore_samples(self.hops_data)

        if len(self.hops_data) == 0:
            raise UserWarning

        log.info("Found {} samples".format(len(self.hops_data)))

        # This type of data isn't 'summarise-able' for general stats, so 
        # skipping straight to heatmap. We also won't write data file to the
        # multiqc_data directory because it would be exactly same as input JSON. 
        self.hops_heatmap()

    def parseJSON(self, f):
        """ Parse the JSON output from HOPS and save the summary statistics """

        try:
            parsed_json = json.load(f['f'])
        except JSONDecodeError as e:
            log.debug("Could not parse HOPS JSON: '{}'".format(f['fn']))
            log.debug(e)
            return None
        
        # Convert JSON to dict for easier manipulation
        for s in parsed_json: 
            s_name = self.clean_s_name(s, f['root'])
            if s_name in self.hops_data: 
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name=s_name)
            self.hops_data[s_name] = {}
            for t in parsed_json[s]:
                self.hops_data[s_name][t] = parsed_json[s][t]
    
    def hops_heatmap(self):
        """ Heatmap showing all statuses for every sample """
        heatmap_numbers = {
            'none': 1,
            'edit_only': 2,
            'damage_only': 3,
            'edit_and_damage': 4
        }

        samples = []
        for s in self.hops_data:
            samples.append(s)

        # As all samples always have same taxa, will take from the first sample
        taxa = []
        for t in self.hops_data[samples[0]]:
            taxa.append(t)

        # Get values from named list into a list of lists required for heatmap
        levels = []
        for s in samples:
            levels.append(self.hops_data[s].values())

        pconfig = {
            'id': 'hops-heatmap',
            'title': 'HOPS: Potential Candidates',
            'xTitle': 'Node',
            'yTitle': 'Sample',
            'min': 0,
            'max': 1,
            'square': False,
            'colstops': [
                [1, '#ffffff'],
                [2, '#ffff33'],
                [3, '#ff7f00'],
                [4, '#e41a1c'],
            ],
            'decimalPlaces': 0,
            'legend': False,
            'datalabels': False,
            'xcats_samples': False,
        }

        self.add_section (
            name = 'Potential Candidates',
            anchor = 'hops_heatmap',
            description = '''
            Heatmap of candidate taxa for downstream aDNA analysis, with 
            intensity representing additive categories of possible 'positive' 
            hits. Note that large numbers of samples can result in Y-axis label 
            overlap, drag down to view all. 

            * Yellow: Edit Distance. 
            * Orange: Damage. 
            * Red: Edit Distance and Damage.
            ''',
            helptext = '''
            HOPS assigns a category based on how many ancient DNA 
            characteristics a given node (i.e. taxon) in a sample has. This goes
            from grey (no characteristics detected), yellow (small edit distance
            from reference), orange (has typical aDNA damage pattern), to
            red (both small edit distance and damage pattern). A red category
            will typically indicates a good candidate for further investigation
            in downstream analysis.

            If data includes many samples, expand plot for full sample list on 
            Y-axis.
            ''',
            plot = heatmap.plot(levels, xcats = taxa, ycats = samples, pconfig = pconfig)
        )
