#!/usr/bin/env python

""" MultiQC module to parse output from QUAST """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='QUAST', anchor='quast',
        href="http://quast.bioinf.spbau.ru/",
        info="is a quality assessment tool for genome assemblies, written by " \
             "the Center for Algorithmic Biotechnology.")

        # Find and load any QUAST reports
        self.quast_data = dict()
        for f in self.find_log_files(config.sp['quast']):
            self.parse_quast_log(f)

        if len(self.quast_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.quast_data)))

        # Write parsed report data to a file
        self.write_data_file(self.quast_data, 'multiqc_quast')

        # Basic Stats Table
        self.quast_general_stats_table()
        
        # Quast Stats Table
        # Only one section, so add to the intro
        self.intro += self.quast_table()
    
    
    def parse_quast_log(self, f):
        lines = f['f'].splitlines()
        
        # Pull out the sample names from the first row
        s_names = lines[0].split("\t")
        for s_name in s_names[1:]:
            if s_name in self.quast_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.quast_data[s_name] = dict()
        
        # Parse remaining stats for each sample
        for l in lines[1:]:
            s = l.split("\t")
            k = s[0]
            for i, v in enumerate(s[1:]):
                s_name = s_names[i+1]
                partials = re.search("(\d+) \+ (\d+) part", v)
                if partials:
                    whole = partials.group(1)
                    partial = partials.group(2)
                    try:
                        self.quast_data[s_name][k] = float(whole)
                        self.quast_data[s_name]["{}_partial".format(k)] = float(partial)
                    except ValueError:
                        self.quast_data[s_name][k] = whole
                        self.quast_data[s_name]["{}_partial".format(k)] = partial
                else:
                    try:
                        self.quast_data[s_name][k] = float(v)
                    except ValueError:
                        self.quast_data[s_name][k] = v
    
    def quast_general_stats_table(self):
        """ Take the parsed stats from the QUAST report and add some to the
        General Statistics table at the top of the report """
        
        headers = OrderedDict()
        headers['N50'] = {
            'title': 'N50 (Kbp)',
            'description': 'N50 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly (kilo base pairs)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.1f}',
            'modify': lambda x: x / 1000
        }
        headers['Total length'] = {
            'title': 'Length (Mbp)',
            'description': 'The total number of bases in the assembly (mega base pairs).',
            'min': 0,
            'suffix': 'bp',
            'scale': 'YlGn',
            'format': '{:.1f}',
            'modify': lambda x: x / 1000000
        }
        self.general_stats_addcols(self.quast_data, headers)
    
    def quast_table(self):
        """ Write some more statistics about the assemblies in a table. """
        
        headers = OrderedDict()
        headers['# misassemblies'] = {
            'title': 'Misassemblies',
            'description': 'The number of positions in the assembled contigs where the left flanking sequence aligns over 1 kbp away from the right flanking sequence on the reference (<i>relocation</i>) or they overlap on more than 1 kbp (<i>relocation</i>) or flanking sequences align on different strands (<i>inversion</i>) or different chromosomes (<i>translocation</i>).</span>',
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:.0f}'
        }
        headers['# mismatches per 100 kbp'] = {
            'title': 'Mismatches',
            'description': 'The number of mismatches per 100 kbp',
            'min': 0,
            'scale': 'YlOrRd',
            'format': '{:.2f}',
        }
        headers['# indels per 100 kbp'] = {
            'title': 'Indels',
            'description': 'The number of indels per 100 kbp',
            'min': 0,
            'scale': 'YlOrRd',
            'format': '{:.2f}',
        }
        headers['# genes'] = {
            'title': 'Genes',
            'description': '# Genes',
            'min': 0,
            'scale': 'YlGnBu',
            'format': '{:.0f}',
            'shared_key': 'gene_count'
        }
        headers['# genes_partial'] = {
            'title': 'Genes (Partial)',
            'description': '# Genes (Partial)',
            'min': 0,
            'scale': 'YlGnBu',
            'format': '{:.0f}',
            'shared_key': 'gene_count'
        }
        headers['Genome fraction (%)'] = {
            'title': 'Genome Fraction',
            'description': 'The total number of aligned bases in the reference, divided by the genome size.',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        return plots.table.plot(self.quast_data, headers)

