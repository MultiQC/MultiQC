#!/usr/bin/env python

""" MultiQC module to parse output from QUAST """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import table, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

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
        for f in self.find_log_files('quast'):
            self.parse_quast_log(f)

        # Filter to strip out ignored sample names
        self.quast_data = self.ignore_samples(self.quast_data)

        if len(self.quast_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.quast_data)))

        # Write parsed report data to a file
        self.write_data_file(self.quast_data, 'multiqc_quast')

        # Basic Stats Table
        self.quast_general_stats_table()

        # Quast Stats Table
        self.add_section (
            name = 'Assembly Statistics',
            anchor = 'quast-stats',
            plot = self.quast_table()
        )
        # Number of contigs plot
        self.add_section (
            name = 'Number of Contigs',
            anchor = 'quast-contigs',
            description = """<p>This plot shows the number of contigs found for each assembly, broken
                    down by length.</p> """,
            plot = self.quast_contigs_barplot()
        )
        # Number of genes plot
        ng_pdata = self.quast_predicted_genes_barplot()
        if ng_pdata:
            self.add_section (
                name = 'Number of Predicted Genes',
                anchor = 'quast-genes',
                description = """<p>This plot shows the number of predicted genes found for each
                          assembly, broken down by length.</p>""",
                plot = ng_pdata
            )

    def parse_quast_log(self, f):
        lines = f['f'].splitlines()

        # Pull out the sample names from the first row
        s_names = lines[0].split("\t")
        # Prepend directory name(s) to sample names as configured
        s_names = [self.clean_s_name(s_name, f['root'])
                   for s_name in s_names]
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
            'modify': lambda x: x / 1000
        }
        headers['Total length'] = {
            'title': 'Length (Mbp)',
            'description': 'The total number of bases in the assembly (mega base pairs).',
            'min': 0,
            'suffix': 'bp',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000
        }
        self.general_stats_addcols(self.quast_data, headers)

    def quast_table(self):
        """ Write some more statistics about the assemblies in a table. """

        headers = OrderedDict()
        headers['N50'] = {
            'title': 'N50 (Kbp)',
            'description': 'N50 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly (kilo base pairs)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'modify': lambda x: x / 1000
        }

        headers['N75'] = {
            'title': 'N75 (Kbp)',
            'description': 'N75 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly (kilo base pairs)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'modify': lambda x: x / 1000
        }
        headers['L50'] = {
            'title': 'L50 (k)',
            'description': 'L50 is the number of contigs larger than N50, i.e. the minimum number of contigs comprising half the total assembly length.',
            'min': 0,
            'suffix': '',
            'scale': 'GnYlRd',
            'modify': lambda x: x / 1000
        }
        headers['L75'] = {
            'title': 'L75 (k)',
            'description': 'L75 is the number of contigs larger than N75, i.e. the minimum number of contigs comprising three quarters of the total assembly length.',
            'min': 0,
            'suffix': '',
            'scale': 'GnYlRd',
            'modify': lambda x: x / 1000
        }
        headers['Largest contig'] = {
            'title': 'Largest contig (Kbp)',
            'description': 'The total number of bases in the assembly (mega base pairs).',
            'min': 0,
            'suffix': 'bp',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000
        }
        headers['Total length'] = {
            'title': 'Length (Mbp)',
            'description': 'The total number of bases in the assembly (mega base pairs).',
            'min': 0,
            'suffix': 'bp',
            'scale': 'YlGn',
            'modify': lambda x: x / 1000000
        }
        headers['# misassemblies'] = {
            'title': 'Misassemblies',
            'description': 'The number of positions in the assembled contigs where the left flanking sequence aligns over 1 kbp away from the right flanking sequence on the reference (relocation) or they overlap on more than 1 kbp (relocation) or flanking sequences align on different strands (inversion) or different chromosomes (translocation).',
            'scale': 'RdYlGn-rev',
            'format': '{,:.0f}'
        }
        headers['# mismatches per 100 kbp'] = {
            'title': 'Mismatches/100kbp',
            'description': 'The number of mismatches per 100 kbp',
            'scale': 'YlOrRd',
            'format': '{:,.2f}',
        }
        headers['# indels per 100 kbp'] = {
            'title': 'Indels/100kbp',
            'description': 'The number of indels per 100 kbp',
            'scale': 'YlOrRd',
            'format': '{:,.2f}',
        }
        headers['# genes'] = {
            'title': 'Genes',
            'description': '# Genes',
            'scale': 'YlGnBu',
            'format': '{:,.0f}',
            'shared_key': 'gene_count'
        }
        headers['# genes_partial'] = {
            'title': 'Genes (Partial)',
            'description': '# Genes (Partial)',
            'scale': 'YlGnBu',
            'format': '{:,.0f}',
            'shared_key': 'gene_count'
        }
        headers['# predicted genes (unique)'] = {
            'title': 'Genes',
            'description': '# Predicted Genes (Unique)',
            'scale': 'YlGnBu',
            'format': '{:,.0f}',
            'shared_key': 'gene_count'
        }
        headers['Genome fraction (%)'] = {
            'title': 'Genome Fraction',
            'description': 'The total number of aligned bases in the reference, divided by the genome size.',
            'max': 100,
            'suffix': '%',
            'scale': 'YlGn'
        }
        config = {
            'id': 'quast_table',
            'namespace': 'QUAST',
            'min': 0,
        }
        return table.plot(self.quast_data, headers, config)

    def quast_contigs_barplot(self):
        """ Make a bar plot showing the number and length of contigs for each assembly """

        # Prep the data
        data = dict()
        for s_name, d in self.quast_data.items():
            try:
                p = dict()
                p['>= 50000 bp'] = d['# contigs (>= 50000 bp)']
                p['25000-50000 bp'] = d['# contigs (>= 25000 bp)'] - d['# contigs (>= 50000 bp)']
                p['10000-25000 bp'] = d['# contigs (>= 10000 bp)'] - d['# contigs (>= 25000 bp)']
                p['5000-10000 bp'] = d['# contigs (>= 5000 bp)'] - d['# contigs (>= 10000 bp)']
                p['1000-5000 bp'] = d['# contigs (>= 1000 bp)'] - d['# contigs (>= 5000 bp)']
                p['0-1000 bp'] = d['# contigs (>= 0 bp)'] - d['# contigs (>= 1000 bp)']
                assert sum(p.values()) == d['# contigs (>= 0 bp)']
                data[s_name] = p
            except AssertionError:
                log.warning("Contig counts didn't add up properly for {}".format(s_name))

        # Define the order of the categories, small to big
        cats = [
            '0-1000 bp',
            '1000-5000 bp',
            '5000-10000 bp',
            '10000-25000 bp',
            '25000-50000 bp',
            '>= 50000 bp',
        ]

        pconfig = {
            'id': 'quast_num_contigs',
            'title': 'Number of Contigs',
            'ylab': '# Contigs',
            'yDecimals': False
        }

        return bargraph.plot(data, cats, pconfig)

    def quast_predicted_genes_barplot(self):
        """
        Make a bar plot showing the number and length of predicted genes
        for each assembly
        """

        # Prep the data
        # extract the ranges given to quast with "--gene-thresholds"
        prefix = '# predicted genes (>= '
        suffix = ' bp)'
        all_thresholds = sorted(list(set([
            int(key[len(prefix):-len(suffix)])
            for _, d in self.quast_data.items()
            for key in d.keys()
            if key.startswith(prefix)
            ])))

        data = {}
        ourpat = '>= {}{} bp'
        theirpat = prefix+"{}"+suffix
        for s_name, d in self.quast_data.items():
            thresholds = sorted(list(set([
                int(key[len(prefix):-len(suffix)])
                for _, x in self.quast_data.items()
                for key in x.keys()
                if key.startswith(prefix)
            ])))
            if len(thresholds)<2: continue

            p = dict()
            try:
                p = { ourpat.format(thresholds[-1],""): d[theirpat.format(thresholds[-1])] }
                for low,high in zip(thresholds[:-1], thresholds[1:]):
                    p[ourpat.format(low,-high)] = d[theirpat.format(low)] - d[theirpat.format(high)]

                assert sum(p.values()) == d[theirpat.format(0)]
            except AssertionError:
                log.warning("Predicted gene counts didn't add up properly for \"{}\"".format(s_name))
            except KeyError:
                log.warning("Not all predicted gene thresholds available for \"{}\"".format(s_name))
            data[s_name] = p

        cats = [ ourpat.format(low,-high if high else "")
                 for low,high in zip(all_thresholds, all_thresholds[1:]+[None]) ]

        if len(cats) > 0:
            return bargraph.plot(data, cats)
        else:
            return None
