#!/usr/bin/env python

""" MultiQC module to parse output from WhatsHap """

from collections import defaultdict, OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    WhatsHap module class, parses WhatsHap output.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
                name='WhatsHap', anchor='whatshap',
                href='https://whatshap.readthedocs.io/',
                info=("""
                 is a program for phasing genomic variants using DNA sequencing
                reads, also called read-based phasing or haplotype assembly. It
                is especially suitable for long reads, but also works well with
                short reads.
                """)
        )

        # Store the whatshap stats results
        self.whatshap_stats = dict()

        # Iterate over all files we found
        for f in self.find_log_files('whatshap/stats', filehandles=True):
            sample, data = self.parse_whatshap_stats(f)
            self.whatshap_stats[sample] = data

        # Filter to strip out ignored sample names
        self.whatshap_stats = self.ignore_samples(self.whatshap_stats)

        # Raise UserWarning if we didn't find any data
        if not self.whatshap_stats:
            raise UserWarning

        # Write parsed report data to a file
        self.write_data_file(self.whatshap_stats, 'multiqc_whatshap_stats')

        # Add whatshap stats to general statistics
        self.whatshap_add_general_stats()

        # Create bargraph with total phased bp
        self.add_bargraph_total_phased()

    def parse_whatshap_stats(self, logfile):
        """ Parse WhatsHap stats file """

        def process_data(data):
            """
            Process the data to convert it to the proper format

            In practice, this only involves the following changes:
                - Convert all numeric values to integers
                - Replace 'nan' values with zero
            """
            for key, value in data.items():
                # Try to make the value an int
                try:
                    data[key] = int(data[key])
                except ValueError:
                    pass

                # Replace 'nan' with zero
                if value == 'nan':
                    data[key] = 0
                
        file_content = logfile['f']
        # This will later be replaced by the sample name from the file, unless
        # we are parsing an empty file
        sample = self.clean_s_name(logfile['s_name'], logfile['root'])
        # Get the header and remove the "#" character from the first field
        header = next(file_content).strip().split()
        header[0] = header[0][1:]

        # The results from this log file, a dictionary for every chromosome
        results = defaultdict(dict)
        # Parse the lines that have the data
        for line in file_content:
            spline = line.strip().split()
            data = {field:value for field, value in zip(header, spline)}
            process_data(data)
            # Remove the sample and chromsome from the data
            sample = data.pop('sample')
            chromosome = data.pop('chromosome')
            # Insert the current line under chromosome
            results[chromosome] = data

        # If we were parsing an empty file, we return dummy data (basically all
        # zero) so the sample statistics are consistent, which eases parsing
        # later on. WhatsHap uses the fictional 'ALL' chromsome to store the
        # summary for all results, so we re-use that to store all zero's
        if not results:
            results = dict()
            results['ALL'] = defaultdict(int)

        return sample, results

    @staticmethod
    def get_summary_field(sample_stats):
        """
        Return the summary field for sample_stats

        If there is only a single chromosome, we use that as the summary field.
        If there are multiple chromosomes, whatshap adds the 'ALL' field which
        contains the summary
        """
        # If there is only a single chromosome we use that as the summary
        # filed. Otherwise, we use the 'ALL' chromosome (it only gets added
        # by WhatsHap stats when there are multiple chromosomes).
        if len(sample_stats) == 1:
            summary_field = list(sample_stats)[0]
        else:
            summary_field = 'ALL'
        return summary_field

    def whatshap_add_general_stats(self):
        """ Add WhatsHap stats to the general statistics table """

        # Store the configuration for each column, this is also used to
        # determine which columns to add to the general statistics table
        # https://whatshap.readthedocs.io/en/latest/guide.html#the-tsv-statistics-format
        general_stats_headers = OrderedDict([
                ('variants', {
                    'id': 'variants',
                    'title': 'Input Variants',
                    'description':
                        'Number of biallelic variants in the input VCF, but '
                        'excluding any non-SNV variants if --only-snvs was '
                        'used',
                    'hidden': False}

                ),
                ('heterozygous_variants', {
                    'id': 'heterozygous_variants',
                    'title': 'Heterozygous Variants',
                    'description':
                        'The number of biallelic, heterozygous variants in '
                        'the input VCF. This is a subset of Input Variants.',
                    'hidden': False}

                ),
                ('heterozygous_snvs', {
                    'id': 'heterozygous_snvs',
                    'title': 'Heterozygous SNVs',
                    'description':
                        'The number of biallelic, heterozygous SNVs in the '
                        'input VCF. This is a subset of Heterozygous '
                        'Variants.',
                    'hidden': False}

                ),
                ('unphased', {
                    'id': 'unphased',
                    'title': 'Unphased Variants',
                    'description':
                        'The number of biallelic, heterozygous variants that '
                        'are not marked as phased in the input VCF. This '
                        'is a subset of heterozygous_variants.',
                    'hidden': False}
                ),
                ('phased', {
                    'id': 'phased',
                    'title': 'Phased Variants',
                    'description':
                        'The number of biallelic, heterozygous variants that '
                        'are marked as phased in the input VCF. This is '
                        'a subset of heterozygous_variants. Also, phased + '
                        'unphased = heterozygous_variants.',
                    'hidden': False}
                ),
                ('phased_snvs', {
                    'id': 'phased_snvs',
                    'title': 'Phased SNVs',
                    'description':
                        'The number of biallelic, heterozygous SNVs that are '
                        'marked as phased in the input VCF. This is a subset '
                        'of phased.',
                    'hidden': False}
                ),
                ('blocks', {
                    'id': 'blocks',
                    'title': 'Blocks',
                    'description':
                        'The total number of phase sets/blocks.',
                    'hidden': False}
                ),
                ('singletons', {
                    'id': 'singletons',
                    'title': 'Singletons',
                    'description':
                        'The number of blocks that contain exactly one '
                        'variant.',
                    'hidden': False}
                ),
                ('bp_per_block_sum', {
                    'id': 'bp_per_block_sum',
                    'title': 'Total Phased bp',
                    'description':
                        'Description of the distribution of non-singleton '
                        'block lengths, where the length of a block is the '
                        'number of basepairs it covers minus 1. That is, a '
                        'block with two variants at positions 2 and 5 has '
                        'length 3.',
                    'hidden': False}
                ),
                ('variant_per_block_avg', {
                    'id': 'variant_per_block_avg',
                    'title': 'Avg Variants per Block',
                    'description':
                        'Description of the distribution of non-singleton '
                        'block sizes, where the size of a block is the number '
                        'of variants it contains. Number of biallelic '
                        'variants in the input VCF, but excluding any non-SNV '
                        'variants if --only-snvs was used.',
                    'hidden': False}
                ),
                ('bp_per_block_avg', {
                    'id': 'bp_per_block_avg',
                    'title': 'Avg bp per Block',
                    'description':
                        'Description of the distribution of non-singleton '
                        'block lengths, where the length of a block is the '
                        'number of basepairs it covers minus 1. That is, a '
                        'block with two variants at positions 2 and 5 has '
                        'length 3.',
                    'hidden': False}
                ),
                ('block_n50', {
                    'id': 'block_n50',
                    'title': 'NG50',
                    'description':
                        'The NG50 value of the distribution of the block '
                        'lengths. Interleaved blocks are cut in order to '
                        'avoid artificially inflating this value.',
                    'hidden': False}
                )
            ])

        general = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)

            # The summary data we are adding to the general statistics
            summary_data = sample_stats[summary_field]

            general_stats = {field: summary_data[field] for field in
                             general_stats_headers}

            general[sample] = general_stats

        self.general_stats_addcols(general, general_stats_headers)

    def add_bargraph_total_phased(self):
        """ Add a bargraph of the total number of phased basepairs """
        pdata = dict()
        for sample, sample_stats in self.whatshap_stats.items():
            # Get the summary field
            summary_field = self.get_summary_field(sample_stats)
            pdata[sample] = { 'Phased basepairs': sample_stats[summary_field]['bp_per_block_sum'] }

        configuration = {
            'id': 'multiqc_whatshap_phased_bp',
            'title': 'Phased basepairs per Sample',
            'anchor': 'multiqc_whatshap_phased_bp'
        }

        self.add_section(
                name='Phased basepairs per Sample',
                anchor='multiqc_whatshap_phased_bp',
                description=
                """
                    This plot show the total number of phased basepairs for
                    each sample.
                """,
                plot = bargraph.plot(pdata))#, configuration))
