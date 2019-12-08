#!/usr/bin/env python

""" MultiQC module to parse output from rna_seqc """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, beeswarm, linegraph, heatmap
from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name="RNA-SeQC", anchor="rna_seqc",
                                            href="https://github.com/broadinstitute/rnaseqc",
                                            info="Fast, efficient RNA-Seq metrics for quality control and process optimization")

        # Check if metrics is version 1 or 2.
        self.rna_seqc_metrics = dict()
        self.rna_seqc2_metrics = dict()
        for file in self.find_log_files('rna_seqc/metrics'):
            self.parse_metrics_rnaseqc(file)

        # detect if using version 2 data
        if len(self.rna_seqc2_metrics) != 0:
            # Filters to strip out ignored sample names
            self.rna_seqc2_metrics = self.ignore_samples(self.rna_seqc2_metrics)

            num_found = len(self.rna_seqc2_metrics)
            log.info("Found {} samples".format(num_found))

            if num_found == 0:
                raise UserWarning

            self.write_data_file(self.rna_seqc2_metrics, 'multiqc_rna_seqc')

            self.rnaseqc2_general_stats()
            self.transcript_associated_plot_v2()
            self.strand_barplot_v2()
            self.bam_statplot()

        if len(self.rna_seqc_metrics) != 0:
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

            num_found = max(len(self.rna_seqc_metrics), len(self.rna_seqc_norm_high_cov),
                            len(self.rna_seqc_norm_medium_cov), len(self.rna_seqc_norm_low_cov))
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

# RNASeQC2 section
    def transcript_associated_plot_v2(self):
        """ Plot a bargraph showing the Transcript-associated reads  """

        # Plot bar graph of groups
        keys = OrderedDict()
        keys['Exonic Rate'] = {'name': 'Exonic', 'color': '#2f7ed8'}
        keys['Intronic Rate'] = {'name': 'Intronic', 'color': '#8bbc21'}
        keys['Intergenic Rate'] = {'name': 'Intergenic', 'color': '#0d233a'}

        # Config for the plot
        pconfig = {
            'id': 'transcript_plot_v2',
            'title': 'RNA-SeQC: Transcript-associated reads',
            'ylab': 'Ratio of Reads',
            'cpswitch': False,
            'ymax': 1,
            'ymin': 0,
            'tt_decimals': 3,
            'cpswitch_c_active': False
        }
        self.add_section(
            name='Transcript associated reads',
            anchor='Transcript_associated_rnaseqc2',
            helptext='All of the above rates are per mapped read. Exonic Rate is the fraction mapping within exons. '
                     'Intronic Rate is the fraction mapping within introns. '
                     'Intergenic Rate is the fraction mapping in the genomic space between genes. ',
            plot=bargraph.plot(self.rna_seqc2_metrics, keys, pconfig)
        )

    def bam_statplot(self):
        pconfig = {
            'id': 'rnaseqc_bam_stat'
        }
        keys = OrderedDict()
        defaults = {
            'min': 0,
            'shared_key': 'read_count',
            'decimalPlaces': 2,
            'modify': lambda x: float(x) / 1000000.0,
        }
        keys['Total Read Number'] = dict(defaults, **{'title': 'Total Read Number'})
        keys['Alternative Alignments'] = dict(defaults, **{'title': 'Alternative Alignments'})
        keys['Chimeric Reads'] = dict(defaults, **{'title': 'Chimeric Reads'})
        keys['Duplicate Reads'] = dict(defaults, **{'title': 'Duplicate Reads'})
        keys['End 1 Mapped Reads'] = dict(defaults, **{'title': 'End 1 Mapped Reads'})
        keys['End 2 Mapped Reads'] = dict(defaults, **{'title': 'End 2 Mapped Reads'})
        keys['End 1 Mismatches'] = dict(defaults, **{'title': 'End 1 Mismatches'})
        keys['End 2 Mismatches'] = dict(defaults, **{'title': 'End 2 Mismatches'})
        keys['End 1 Sense'] = dict(defaults, **{'title': 'End 1 Sense'})
        keys['End 2 Sense'] = dict(defaults, **{'title': 'End 2 Sense'})
        keys['Ambiguous Reads'] = dict(defaults, **{'title': 'Ambiguous Reads'})
        keys['High Quality Reads'] = dict(defaults, **{'title': 'High Quality Reads'})
        keys['Low Quality Reads'] = dict(defaults, **{'title': 'Low Quality Reads'})
        keys['Mapped Duplicate Reads'] = dict(defaults, **{'title': 'Mapped Duplicate Reads'})
        keys['Mapped Reads'] = dict(defaults, **{'title': 'Mapped Reads'})
        keys['Mapped Unique Reads'] = dict(defaults, **{'title': 'Mapped Unique Reads'})
        keys['Non-Globin Reads'] = dict(defaults, **{'title': 'Non-Globin Reads'})
        keys['Non-Globin Duplicate Reads'] = dict(defaults, **{'title': 'Non-Globin Duplicate Reads'})
        keys['rRNA Reads'] = dict(defaults, **{'title': 'rRNA Reads'})
        keys['Unique Mapping, Vendor QC Passed Reads'] = dict(defaults,
                                                              **{'title': 'Unique Mapping, Vendor QC Passed Reads'})

        self.add_section(
            name='Read Counts',
            anchor='rnaseqc-bam_stat',
            description='All numbers are reported in millions.',
            plot=beeswarm.plot(self.rna_seqc2_metrics, keys, pconfig)
        )

    def strand_barplot_v2(self):
        """ Plot a bargraph showing the strandedness of alignments """
        # Plot bar graph of groups
        keys = ['End 1 Sense', 'End 1 Antisense', 'End 2 Sense', 'End 2 Antisense']
        # Config for the plot
        pconfig = {
            'id': 'strand_spec_barplot',
            'title': 'RNA-SeQC: Strand Specificity',
            'ylab': '% Reads',
            'cpswitch_counts_label': '# Reads',
            'cpswitch_percent_label': '% Reads',
            'ymin': 0,
            'cpswitch_c_active': False
        }
        self.add_section(
            name='Strand Specificity',
            anchor='rna_seqc_strand_specificity_rnaseqc2',
            helptext='End 1/2 Sense are the number of End 1 or 2 reads that were sequenced in the sense direction. '
                     'Similarly, End 1/2 Antisense are the number of End 1 or 2 reads that were sequenced in the '
                     'antisense direction',
            plot=bargraph.plot(self.rna_seqc2_metrics, keys, pconfig)
        )

    def rnaseqc2_general_stats(self):
        """
        Add alignment rate to the general stats table
        """
        headers = OrderedDict()
        headers['Expression Profiling Efficiency'] = {
            'namespace': 'RNA-SeQCv2',
            'title': '% Expression Efficiency',
            'description': 'Expression Profiling Efficiency: Ratio of exon reads to total reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0
        }
        headers['Genes Detected'] = {
            'namespace': 'RNA-SeQCv2',
            'title': '# Genes',
            'description': 'Number of genes detected with at least 5 reads.',
            'min': 0,
            'scale': 'Bu',
            'format': '{:,.0f}'
        }
        headers['rRNA Rate'] = {
            'namespace': 'RNA-SeQCv2',
            'title': '% rRNA Alignment',
            'description': ' rRNA reads (non-duplicate and duplicate reads) per total reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'modify': lambda x: float(x) * 100.0
        }

        self.general_stats_addcols(self.rna_seqc2_metrics, headers)

    def parse_metrics_rnaseqc(self, f):
        """
        Parse the metrics.tsv file from RNA-SeQC
        """
        headers = list()

        lines = f['f'].splitlines()
        s = lines[0].split('\t')
        version = 0
        attempt = None
        try:
            attempt = s[2]
            version = 1
        except IndexError:
            version = 2

        if version == 1:
            headers = None
            for l in f['f'].splitlines():
                s = l.strip().split("\t")
                if headers is None:
                    headers = s
                else:
                    s_name = s[headers.index('Sample')]
                    data = dict()
                    for idx, h in enumerate(headers):
                        try:
                            data[h] = float(s[idx])
                        except ValueError:
                            data[h] = s[idx]
                    self.rna_seqc_metrics[s_name] = data
        elif version == 2:
            # handle header creation (get first column from the file).
            for l in f['f'].splitlines():
                s = l.split('\t')
                if s[0] == "Total Reads":
                    s[0] = "Total Read Number"
                headers.append(s[0])

            # sample name
            s_name = f['f'].splitlines()[0].split('\t')[1].split('.bam')[0]

            data = dict()
            i = 0
            for l in f['f'].splitlines():
                if i == 0:
                    data[headers[i]] = s_name
                else:
                    s = l.split('\t')
                    data[headers[i]] = s[1]
                i += 1

            self.rna_seqc2_metrics[s_name] = data



    def rnaseqc_general_stats(self):
        """
        Add alignment rate to the general stats table
        """
        headers = OrderedDict()
        headers['Expression Profiling Efficiency'] = {
            'namespace': 'RNA-SeQCv1',
            'title': '% Expression Efficiency',
            'description': 'Expression Profiling Efficiency: Ratio of exon reads to total reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0
        }
        headers['Genes Detected'] = {
            'namespace': 'RNA-SeQCv1',
            'title': '# Genes',
            'description': 'Number of genes detected with at least 5 reads.',
            'min': 0,
            'scale': 'Bu',
            'format': '{:,.0f}'
        }
        headers['rRNA rate'] = {
            'namespace': 'RNA-SeQCv1',
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
            'id': 'transcript_plot_v1',
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
            anchor = 'Transcript_associated_rnaseqc1',
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
            'id': 'strand_spec_v1',
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
        data = list()
        data_labels = list()
        if len(self.rna_seqc_norm_high_cov) > 0:
            data.append(self.rna_seqc_norm_high_cov)
            data_labels.append({'name': 'High Expressed'})
        if len(self.rna_seqc_norm_medium_cov) > 0:
            data.append(self.rna_seqc_norm_medium_cov)
            data_labels.append({'name': 'Medium Expressed'})
        if len(self.rna_seqc_norm_low_cov) > 0:
            data.append(self.rna_seqc_norm_low_cov)
            data_labels.append({'name': 'Low Expressed'})
        pconfig = {
            'id': 'gene_body_coverage_rnaseqc',
            'title': 'RNA-SeQC: Gene Body Coverage',
            'ylab': '% Coverage',
            'xlab': "Gene Body Percentile (5' -> 3')",
            'xmin': 0,
            'xmax': 100,
            'tt_label': "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
            'data_labels': data_labels
        }
        if len(data) > 0:
            self.add_section (
                name = 'Gene Body Coverage',
                anchor = 'rseqc-rna_seqc_mean_coverage',
                helptext = 'The metrics are calculated across the transcripts with tiered expression levels.',
                plot = linegraph.plot(data, pconfig)
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
                'id': 'plot_correlation_rnaseqc',
                'title': 'RNA-SeQC: {} Sample Correlation'.format(corr_type)
            }
            self.add_section (
                name = '{} Correlation'.format(corr_type),
                anchor = 'rseqc-rna_seqc_correlation',
                plot = heatmap.plot(data[1], data[0], data[0], pconfig)
            )
