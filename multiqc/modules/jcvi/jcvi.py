#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='JCVI Genome Annotation', anchor='jcvi',
                                            href="https://pypi.org/project/jcvi/",
                                            info="computes statistics on genome annotation.")

        # Parse logs
        self.jcvi = dict()
        for f in self.find_log_files('jcvi', filehandles=True):
            self.parse_jcvi(f)

        # Filter to strip out ignored sample names
        self.jcvi = self.ignore_samples(self.jcvi)

        if len(self.jcvi) == 0:
            raise UserWarning

        log.info("Found {} logs".format(len(self.jcvi)))
        self.write_data_file(self.jcvi, 'multiqc_jcvi')

        # Add most important JCVI annotation stats to the general table
        headers = OrderedDict()
        headers['genes'] = {
            'title': 'Number of genes',
            'description': 'Number of genes',
            'format': '{:i}',
        }
        headers['transcripts'] = {
            'title': 'Number of transcripts',
            'description': 'Number of transcripts',
            'format': '{:i}',
        }
        headers['mean_gene_size'] = {
            'title': 'Mean gene length (bp)',
            'description': 'Mean gene length',
            'format': '{:i}',
        }
        self.general_stats_addcols(self.jcvi, headers)

        self.add_section(plot=self.jcvi_barplot_genes(),
                         name="Number of genes",
                         description="Total number of genes found in each dataset")

        self.add_section(plot=self.jcvi_barplot_transcripts(),
                         name="Number of transcripts",
                         description="Total number of transcripts found in each dataset")

        self.add_section(plot=self.jcvi_barplot_exons(),
                         name="Number of exons",
                         description="Total number of exons found in each dataset")

        self.add_section(plot=self.jcvi_barplot_genes_len(),
                         name="Mean size of genes",
                         description="Mean size of all genes foudn in each dataset")

        self.add_section(plot=self.jcvi_barplot_transcripts_len(),
                         name="Mean size of transcripts",
                         description="Mean size of all transcripts foudn in each dataset")

        self.add_section(plot=self.jcvi_barplot_exons_len(),
                         name="Mean size of exons",
                         description="Mean size of all exons foudn in each dataset")

        self.add_section(plot=self.jcvi_barplot_transcripts_per_genes(),
                         name="Transcripts per gene",
                         description="Mean and maximum number of transcripts per gene")

        self.add_section(plot=self.jcvi_barplot_isoforms(),
                         name="Isoforms",
                         description="Number of genes found to have multiple isoforms")

        self.add_section(plot=self.jcvi_barplot_exons_per_genes(),
                         name="Exons per gene",
                         description="Mean number of exnos per gene")

        gene_length_plot = self.jcvi_linegraph_gene_length()
        if gene_length_plot:
            self.add_section(
                name='Gene length distribution',
                anchor='jcvi_gene_length',
                description='This plot shows the distribution of gene length.',
                helptext='''
                Values are binned in buckets of 100bp (e.g. the value at 150bp represents the number of genes having a length between 100 and 199bp).
                ''',
                plot=gene_length_plot
            )

        exon_length_plot = self.jcvi_linegraph_exon_length()
        if exon_length_plot:
            self.add_section(
                name='Exon length distribution',
                anchor='jcvi_exon_length',
                description='This plot shows the distribution of exon length.',
                helptext='''
                Values are binned in buckets of 25bp (e.g. the value at 112bp represents the number of exons having a length between 100 and 124bp).
                ''',
                plot=exon_length_plot
            )

        intron_length_plot = self.jcvi_linegraph_intron_length()
        if intron_length_plot:
            self.add_section(
                name='Intron length distribution',
                anchor='jcvi_intron_length',
                description='This plot shows the distribution of intron length.',
                helptext='''
                Values are binned in buckets of 25bp (e.g. the value at 112bp represents the number of introns having a length between 100 and 124bp).
                ''',
                plot=intron_length_plot
            )

        exon_count_plot = self.jcvi_linegraph_exon_count()
        if exon_count_plot:
            self.add_section(
                name='Exon count distribution',
                anchor='jcvi_exon_count',
                description='This plot shows the distribution of exon number per genes.',
                helptext='''
                If you look at y-axis value corresponding to the x-axis value "3", you get the number of genes having 3 exons.
                ''',
                plot=exon_count_plot
            )

    def jcvi_linegraph_gene_length(self):
        plot_config = {
            'id': 'jcvi_gene_length_plot',
            'title': 'JCVI: Gene length repartition',
            'ylab': '# genes',
            'xlab': 'Gene length (bp)',
            'xDecimals': False
        }

        plot_data = {x: self.jcvi[x]['gene_length'] for x in self.jcvi if 'gene_length' in self.jcvi[x]}

        if not plot_data:
            return None

        return linegraph.plot(plot_data, plot_config)

    def jcvi_linegraph_exon_length(self):
        plot_config = {
            'id': 'jcvi_exon_length_plot',
            'title': 'JCVI: Exon length repartition',
            'ylab': '# exons',
            'xlab': 'Exon length (bp)',
            'xDecimals': False
        }

        plot_data = {x: self.jcvi[x]['exon_length'] for x in self.jcvi if 'exon_length' in self.jcvi[x]}

        if not plot_data:
            return None

        return linegraph.plot(plot_data, plot_config)

    def jcvi_linegraph_intron_length(self):
        plot_config = {
            'id': 'jcvi_intron_length_plot',
            'title': 'JCVI: Intron length repartition',
            'ylab': '# introns',
            'xlab': 'Intron length (bp)',
            'xDecimals': False
        }

        plot_data = {x: self.jcvi[x]['intron_length'] for x in self.jcvi if 'intron_length' in self.jcvi[x]}

        if not plot_data:
            return None

        return linegraph.plot(plot_data, plot_config)

    def jcvi_linegraph_exon_count(self):
        plot_config = {
            'id': 'jcvi_exon_count_plot',
            'title': 'JCVI: Exon count repartition',
            'ylab': '# genes',
            'xlab': 'Exon count',
            'xDecimals': False,
            'ymin': 0,
        }

        plot_data = {x: self.jcvi[x]['exon_count'] for x in self.jcvi if 'exon_count' in self.jcvi[x]}

        if not plot_data:
            return None

        return linegraph.plot(plot_data, plot_config)

    def parse_jcvi(self, f):

        s_name = None

        # Look at the first three lines, they are always the same
        first_line = f['f'].readline()
        second_line = f['f'].readline()
        third_line = f['f'].readline()
        # If any of these fail, it's probably not a jcvi summary file
        if not all((first_line.startswith("==================================================================="),
                   second_line.startswith("                                                   o            all"),
                   third_line.startswith("-------------------------------------------------------------------"))):
            return

        s_name = f['s_name']

        if s_name in self.jcvi:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))

        self.jcvi[s_name] = dict()

        match = re.search(r'Max number of transcripts per gene\s+([\d,]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['transcripts_per_gene'] = int(match.group(1).replace(',', ''))

        match = re.search(r'Mean exon size\s+([\d,.]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['mean_exons_size'] = float(match.group(1).replace(',', ''))

        match = re.search(r'Mean gene locus size \(first to last exon\)\s+([\d,.]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['mean_gene_size'] = float(match.group(1).replace(',', ''))

        match = re.search(r'Mean number of distinct exons per gene\s+([\d,.]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['mean_exon_number'] = float(match.group(1).replace(',', ''))

        match = re.search(r'Mean number of transcripts per gene\s+([\d,.]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['mean_transcript_number'] = float(match.group(1).replace(',', ''))

        match = re.search(r'Mean transcript size \(UTR, CDS\)\s+([\d,.]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['mean_transcript_size'] = float(match.group(1).replace(',', ''))

        match = re.search(r'Number of distinct exons\s+([\d,]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['exons'] = int(match.group(1).replace(',', ''))

        match = re.search(r'Number of genes\s+([\d,]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['genes'] = int(match.group(1).replace(',', ''))

        match = re.search(r'Number of genes with alternative transcript variants\s+(\d+)\s\(([\d.]+)%\)', f['f'].readline())
        if match:
            self.jcvi[s_name]['genes_with_alt'] = int(match.group(1).replace(',', ''))
            self.jcvi[s_name]['genes_with_alt_percent'] = float(match.group(2).replace(',', ''))

        match = re.search(r'Number of multi-exon genes\s+([\d,]+)\s\(([\d.]+)%\)', f['f'].readline())
        if match:
            self.jcvi[s_name]['multiexon_genes'] = int(match.group(1).replace(',', ''))
            self.jcvi[s_name]['multiexon_genes_percent'] = float(match.group(2).replace(',', ''))

        match = re.search(r'Number of predicted transcripts\s+([\d,]+)', f['f'].readline())
        if match:
            self.jcvi[s_name]['transcripts'] = int(match.group(1).replace(',', ''))

        match = re.search(r'Number of single-exon genes\s+([\d,]+)\s\(([\d.]+)%\)', f['f'].readline())
        if match:
            self.jcvi[s_name]['singleexon_genes'] = int(match.group(1).replace(',', ''))
            self.jcvi[s_name]['singleexon_genes_percent'] = float(match.group(2).replace(',', ''))

        exon_count = self.parse_hists(os.path.join(f['root'], 'Exon_Count', s_name + '.txt'))
        if exon_count:
            self.jcvi[s_name]['exon_count'] = exon_count

        exon_length = self.parse_hists(os.path.join(f['root'], 'Exon_Length', s_name + '.txt'), bin_by=25)
        if exon_length:
            self.jcvi[s_name]['exon_length'] = exon_length

        gene_length = self.parse_hists(os.path.join(f['root'], 'Gene_Length', s_name + '.txt'), bin_by=100)
        if gene_length:
            self.jcvi[s_name]['gene_length'] = gene_length

        intron_length = self.parse_hists(os.path.join(f['root'], 'Intron_Length', s_name + '.txt'), bin_by=25)
        if intron_length:
            self.jcvi[s_name]['intron_length'] = intron_length

        self.add_data_source(f, s_name)

    def parse_hists(self, stat_file, bin_by=1):

        stat_table = {}

        if os.path.isfile(stat_file):
            with open(stat_file, 'r') as stat_handle:
                for stat_line in stat_handle:
                    val = ((int(stat_line) / bin_by) * bin_by) + bin_by / 2
                    if val not in stat_table:
                        stat_table[val] = 1
                    else:
                        stat_table[val] += 1

        # Add missing values
        if stat_table:
            values = sorted(stat_table.keys())
            for required_val in range(int(bin_by / 2), int(max(values)), bin_by):
                if required_val not in stat_table:
                    stat_table[required_val] = 0

        return stat_table

    def jcvi_barplot_genes(self):
        keys = OrderedDict()
        keys['multiexon_genes'] = {'name': 'Multi-exon'}
        keys['singleexon_genes'] = {'name': 'Single-exon'}

        plot_config = {
            'id': 'jcvi_plot_singlemultiexons',
            'title': 'JCVI: Number of genes',
            'ylab': '# Counts',
            'cpswitch_counts_label': 'Number of genes'
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_transcripts(self):
        keys = OrderedDict()
        keys['transcripts'] = {'name': 'transcripts'}

        plot_config = {
            'id': 'jcvi_plot_transcripts',
            'title': 'JCVI: Number of transcripts',
            'ylab': '# Counts',
            'cpswitch': False
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_exons(self):
        keys = OrderedDict()
        keys['exons'] = {'name': 'exons'}

        plot_config = {
            'id': 'jcvi_plot_exons',
            'title': 'JCVI: Number of exons',
            'ylab': '# Counts',
            'cpswitch': False
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_genes_len(self):
        keys = OrderedDict()
        keys['mean_gene_size'] = {'name': 'mean_gene_size'}

        plot_config = {
            'id': 'jcvi_plot_genes_len',
            'title': 'JCVI: Mean size of genes',
            'ylab': '# Counts',
            'cpswitch': False,
            'tt_decimals': 1
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_transcripts_len(self):
        keys = OrderedDict()
        keys['mean_transcript_size'] = {'name': 'mean_transcript_size'}

        plot_config = {
            'id': 'jcvi_plot_transcripts_len',
            'title': 'JCVI: Mean size of transcripts',
            'ylab': '# Counts',
            'cpswitch': False,
            'tt_decimals': 1
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_exons_len(self):
        keys = OrderedDict()
        keys['mean_exons_size'] = {'name': 'mean_exons_size'}

        plot_config = {
            'id': 'jcvi_plot_exons_len',
            'title': 'JCVI: Mean size of exons',
            'ylab': '# Counts',
            'cpswitch': False,
            'tt_decimals': 1
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_transcripts_per_genes(self):
        keys = OrderedDict()
        keys['mean_transcript_number'] = {'name': 'Mean transcripts per genes'}
        keys['transcripts_per_gene'] = {'name': 'Maximum transcripts per genes'}

        plot_config = {
            'id': 'jcvi_plot_transcripts_per_genes',
            'title': 'JCVI: Transcripts per gene',
            'ylab': '# transcripts per gene',
            'cpswitch': False,
            'tt_decimals': 1
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_isoforms(self):
        keys = OrderedDict()
        keys['genes_with_alt'] = {'name': 'Genes with multiple isoforms'}

        plot_config = {
            'id': 'jcvi_plot_isoforms',
            'title': 'JCVI: Genes with multiple isoforms',
            'ylab': '# genes with multiple isoforms',
            'cpswitch': False
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_barplot_exons_per_genes(self):
        keys = OrderedDict()
        keys['mean_exon_number'] = {'name': 'Mean exons per genes'}

        plot_config = {
            'id': 'jcvi_plot_exons_per_genes',
            'title': 'JCVI: Exons per gene',
            'ylab': '# exons per gene',
            'cpswitch': False,
            'tt_decimals': 1
        }

        return bargraph.plot(self.jcvi, keys, plot_config)
