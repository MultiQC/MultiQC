#!/usr/bin/env python

""" MultiQC module to parse output from RNA-SeQC """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, heatmap

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='RNA-SeQC', anchor='rna_seqc',
        href='http://archive.broadinstitute.org/cancer/cga/rna_seqc',
        info="is a java program which computes a series of quality control metrics for RNA-seq data.")

        # Parse metrics information.
        self.rna_seqc_metrics = dict()
        for f in self.find_log_files('rna_seqc/metrics'):
            self.parse_metrics(f)

        # Parse normalised coverage information.
        self.rna_seqc_norm_high_cov = dict()
        self.rna_seqc_norm_medium_cov = dict()
        self.rna_seqc_norm_low_cov = dict()
        for f in self.find_log_files('rna_seqc/coverage'):
            self.parse_coverage(f)

        # Parse correlation matrices
        self.rna_seqc_pearson = None
        self.rna_seqc_spearman = None
        for f in self.find_log_files('rna_seqc/correlation'):
            self.parse_correlation(f)

        # Filters to strip out ignored sample names
        self.rna_seqc_metrics = self.ignore_samples(self.rna_seqc_metrics)
        self.rna_seqc_norm_high_cov = self.ignore_samples(self.rna_seqc_norm_high_cov)
        self.rna_seqc_norm_medium_cov = self.ignore_samples(self.rna_seqc_norm_medium_cov)
        self.rna_seqc_norm_low_cov = self.ignore_samples(self.rna_seqc_norm_low_cov)
        # TODO: self.rna_seqc_pearson and self.rna_seqc_spearman are trickier to filter

        num_found = max( len(self.rna_seqc_metrics), len(self.rna_seqc_norm_high_cov),
                         len(self.rna_seqc_norm_medium_cov), len(self.rna_seqc_norm_low_cov) )
        if self.rna_seqc_pearson is not None:
            num_found += 1
        if self.rna_seqc_spearman is not None:
            num_found += 1
        if num_found == 0:
            raise UserWarning

        log.info("Found {} samples".format(num_found))
        self.write_data_file(self.rna_seqc_metrics, 'multiqc_rna_seqc')

        self.rnaseqc_general_stats()
        self.transcript_associated_plot()
        self.plot_correlation_heatmap()
        self.strand_barplot()
        self.coverage_lineplot()

    def parse_metrics(self, f):
        """
        Parse the metrics.tsv file from RNA-SeQC
        """
        headers = None
        for l in f['f'].splitlines():
            s = l.strip().split("\t")
            if headers is None:
                headers = s
            else:
                s_name = s[ headers.index('Sample') ]
                data = dict()
                for idx, h in enumerate(headers):
                    try:
                        data[h] = float(s[idx])
                    except ValueError:
                        data[h] = s[idx]
                self.rna_seqc_metrics[s_name] = data

    def rnaseqc_general_stats (self):
        """
        Add alignment rate to the general stats table
        """
        headers = OrderedDict()
        headers['Expression Profiling Efficiency'] = {
            'title': '% Expression Efficiency',
            'description': 'Expression Profiling Efficiency: Ratio of exon reads to total reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0
        }
        headers['Genes Detected'] = {
            'title': '# Genes',
            'description': 'Number of genes detected with at least 5 reads.',
            'min': 0,
            'scale': 'Bu',
            'format': '{:,.0f}'
        }
        headers['rRNA rate'] = {
            'title': '% rRNA Alignment',
            'description': ' rRNA reads (non-duplicate and duplicate reads) per total reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'modify': lambda x: float(x) * 100.0
        }

        self.general_stats_addcols(self.rna_seqc_metrics, headers)

    def transcript_associated_plot (self):
        """ Plot a bargraph showing the Transcript-associated reads  """

        # Plot bar graph of groups
        keys = OrderedDict()
        keys['Exonic Rate'] = { 'name': 'Exonic', 'color': '#2f7ed8' }
        keys['Intronic Rate'] = { 'name': 'Intronic', 'color': '#8bbc21' }
        keys['Intergenic Rate'] = { 'name': 'Intergenic', 'color': '#0d233a'}

        # Config for the plot
        pconfig = {
            'id': 'rna_seqc_position_plot',
            'title': 'RNA-SeQC: Transcript-associated reads',
            'ylab': 'Ratio of Reads',
            'cpswitch': False,
            'ymax': 1,
            'ymin': 0,
            'tt_decimals': 3,
            'cpswitch_c_active': False
        }
        self.add_section (
            name = 'Transcript-associated reads',
            anchor = 'Transcript_associated',
            helptext = 'All of the above rates are per mapped read. Exonic Rate is the fraction mapping within exons. '
                       'Intronic Rate is the fraction mapping within introns. '
                       'Intergenic Rate is the fraction mapping in the genomic space between genes. ',
            plot = bargraph.plot(self.rna_seqc_metrics, keys, pconfig)
        )



    def strand_barplot(self):
        """ Plot a bargraph showing the strandedness of alignments """
        # Plot bar graph of groups
        keys = [ 'End 1 Sense', 'End 1 Antisense', 'End 2 Sense', 'End 2 Antisense' ]
        # Config for the plot
        pconfig = {
            'id': 'rna_seqc_strandedness_plot',
            'title': 'RNA-SeQC: Strand Specificity',
            'ylab': '% Reads',
            'cpswitch_counts_label': '# Reads',
            'cpswitch_percent_label': '% Reads',
            'ymin': 0,
            'cpswitch_c_active': False
        }
        self.add_section (
            name = 'Strand Specificity',
            anchor = 'rna_seqc_strand_specificity',
            helptext = 'End 1/2 Sense are the number of End 1 or 2 reads that were sequenced in the sense direction. '
                       'Similarly, End 1/2 Antisense are the number of End 1 or 2 reads that were sequenced in the '
                       'antisense direction',
            plot = bargraph.plot(self.rna_seqc_metrics, keys, pconfig)
        )


    def parse_coverage (self, f):
        """ Parse the RNA-SeQC Normalised Coverage Files """
        data = dict()
        s_names = None
        j = 1
        for l in f['f'].splitlines():
            s = l.strip().split("\t")
            if s_names is None:
                s_names = s
                for s_name in s_names:
                    data[s_name] = dict()
            else:
                for i, v in enumerate(s):
                    data[s_names[i]][j] = float(v)
                j += 1
        if f['fn'] == 'meanCoverageNorm_high.txt':
            self.rna_seqc_norm_high_cov.update(data)
        elif f['fn'] == 'meanCoverageNorm_medium.txt':
            self.rna_seqc_norm_medium_cov.update(data)
        elif f['fn'] == 'meanCoverageNorm_low.txt':
            self.rna_seqc_norm_low_cov.update(data)

    def coverage_lineplot (self):
        """ Make HTML for coverage line plots """
        # Add line graph to section
        pconfig = {
            'id': 'rna_seqc_mean_coverage_plot',
            'title': 'RNA-SeQC: Gene Body Coverage',
            'ylab': '% Coverage',
            'xlab': "Gene Body Percentile (5' -> 3')",
            'xmin': 0,
            'xmax': 100,
            'tt_label': "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
            'data_labels': [
                {'name': 'High Expressed'},
                {'name': 'Medium Expressed'},
                {'name': 'Low Expressed'}
            ]
        }
        self.add_section (
            name = 'Gene Body Coverage',
            anchor = 'rseqc-rna_seqc_mean_coverage',
            helptext = 'The metrics are calculated across the transcripts with tiered expression levels.',
            plot = linegraph.plot( [
                self.rna_seqc_norm_high_cov,
                self.rna_seqc_norm_medium_cov,
                self.rna_seqc_norm_low_cov
                ], pconfig)
        )

    def parse_correlation(self, f):
        """ Parse RNA-SeQC correlation matrices """
        s_names = None
        data = list()
        for l in f['f'].splitlines():
            s = l.strip().split("\t")
            if s_names is None:
                s_names = [ x for x in s if x != '' ]
            else:
                data.append(s[1:])
        if f['fn'] == 'corrMatrixPearson.txt':
            self.rna_seqc_pearson = (s_names, data)
        elif f['fn'] == 'corrMatrixSpearman.txt':
            self.rna_seqc_spearman = (s_names, data)


    def plot_correlation_heatmap(self):
        """ Return HTML for correlation heatmap """
        data = None
        corr_type = None
        correlation_type = getattr(config, 'rna_seqc' ,{}).get('default_correlation', 'spearman')
        if self.rna_seqc_spearman is not None and correlation_type != 'pearson':
            data = self.rna_seqc_spearman
            corr_type = 'Spearman'
        elif self.rna_seqc_pearson is not None:
            data = self.rna_seqc_pearson
            corr_type = 'Pearson'
        if data is not None:
            pconfig = {
                'id': 'rna_seqc_correlation_heatmap',
                'title': 'RNA-SeQC: {} Sample Correlation'.format(corr_type)
            }
            self.add_section (
                name = '{} Correlation'.format(corr_type),
                anchor = 'rseqc-rna_seqc_correlation',
                plot = heatmap.plot(data[1], data[0], data[0], pconfig)
            )
