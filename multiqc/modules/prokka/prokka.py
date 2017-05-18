#!/usr/bin/env python

""" MultiQC module to parse output from Prokka """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Prokka', anchor='prokka',
        href='http://www.vicbioinformatics.com/software.prokka.shtml',
        info="is a software tool for the rapid annotation of prokaryotic genomes.")

        # Parse logs
        self.prokka = dict()
        for f in self.find_log_files('prokka', filehandles=True):
            self.parse_prokka(f)

        # Filter to strip out ignored sample names
        self.prokka = self.ignore_samples(self.prokka)

        if len(self.prokka) == 0:
            log.debug("Could not find any Prokka data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.prokka)))
        self.write_data_file(self.prokka, 'multiqc_prokka')

        # Add most important Prokka annotation stats to the general table
        headers = OrderedDict()
        headers['organism'] = {
            'title': 'Organism',
            'description': 'Organism',
        }
        headers['contigs'] = {
            'title': 'Contigs',
            'description': 'Number of contigs',
            'min': 0,
        }
        headers['bases'] = {
            'title': 'Bases',
            'description': 'Number of bases',
            'min': 0,
            'format': '{:i}%',
            'hidden': True
        }
        headers['CDS'] = {
            'title': 'CDS',
            'description': 'Number of CDS',
            'min': 0,
            'format': '{:i}%',
        }
        self.general_stats_addcols(self.prokka, headers)

        # User can set configuration attributes, 'prokka_table', and
        # 'prokka_barplot', to specify whether to include a table or a barplot, or both.
        # Default is to make a plot only.
        if getattr(config, 'prokka_table', False):
            self.add_section( plot = self.prokka_table() )
        if getattr(config, 'prokka_barplot', True):
            self.add_section( plot = self.prokka_barplot() )


    def parse_prokka(self, f):
        """ Parse prokka txt summary files.

        Prokka summary files are difficult to identify as there are practically
        no distinct prokka identifiers in the filenames or file contents. This
        parser makes an attempt using the first three lines, expected to contain
        organism, contigs, and bases statistics.
        """

        s_name = None

        # Look at the first three lines, they are always the same
        first_line = f['f'].readline()
        contigs_line = f['f'].readline()
        bases_line = f['f'].readline()
        # If any of these fail, it's probably not a prokka summary file
        if not all((first_line.startswith("organism:"),
                   contigs_line.startswith("contigs:"),
                   bases_line.startswith("bases:"))):
            return

        # Get organism and sample name from the first line
        # Assumes organism name only consists of two words,
        # i.e. 'Genusname speciesname', and that the remaining
        # text on the organism line is the sample name.
        organism = " ".join(first_line.strip().split(":", 1)[1].split()[:2])
        s_name = " ".join(first_line.split()[3:])
        self.prokka[s_name] = dict()
        self.prokka[s_name]['organism'] = organism
        self.prokka[s_name]['contigs'] = int(contigs_line.split(":")[1])
        self.prokka[s_name]['bases'] = int(bases_line.split(":")[1])

        # Get additional info from remaining lines
        for line in f['f']:
            description, value = line.split(":")
            try:
                self.prokka[s_name][description] = int(value)
            except ValueError:
                log.warning("Unable to parse line: '%s'", line)

        self.add_data_source(f, s_name)


    def prokka_table(self):
        """ Make basic table of the annotation stats """

        # Specify the order of the different possible categories
        headers = OrderedDict()
        headers['organism'] = {
                'title': 'Organism',
                'description': 'Organism name',
        }
        headers['contigs'] = {
                'title': '# contigs',
                'description': 'Number of contigs in assembly',
                'format': '{:i}',
        }
        headers['bases'] = {
                'title': '# bases',
                'description': 'Number of nucleotide bases in assembly',
                'format': '{:i}',
        }
        headers['CDS'] = {
                'title': '# CDS',
                'description': 'Number of annotated CDS',
                'format': '{:i}',
        }
        headers['rRNA'] = {
                'title': '# rRNA',
                'description': 'Number of annotated rRNA',
                'format': '{:i}',
        }
        headers['tRNA'] = {
                'title': '# tRNA',
                'description': 'Number of annotated tRNA',
                'format': '{:i}',
        }
        headers['tmRNA'] = {
                'title': '# tmRNA',
                'description': 'Number of annotated tmRNA',
                'format': '{:i}',
        }
        headers['misc_RNA'] = {
                'title': '# misc RNA',
                'description': 'Number of annotated misc. RNA',
                'format': '{:i}',
        }
        headers['sig_peptide'] = {
                'title': '# sig_peptide',
                'description': 'Number of annotated sig_peptide',
                'format': '{:i}',
        }

        table_config = {
            'namespace': 'prokka',
            'min': 0,
        }

        return table.plot(self.prokka, headers, table_config)

    def prokka_barplot(self):
        """ Make a basic plot of the annotation stats """

        # Specify the order of the different categories
        keys = OrderedDict()
        keys['CDS'] =           { 'name': 'CDS' }
        keys['rRNA'] =          { 'name': 'rRNA' }
        keys['tRNA'] =          { 'name': 'tRNA' }
        keys['tmRNA'] =         { 'name': 'tmRNA' }
        keys['misc_RNA'] =      { 'name': 'misc RNA' }
        keys['sig_peptide'] =   { 'name': 'Signal peptides' }

        plot_config = {
            'id': 'prokka_plot',
            'title': 'Prokka',
            'ylab': '# Counts',
            'cpswitch_counts_label': 'Features'
        }

        return bargraph.plot(self.prokka, keys, plot_config)
