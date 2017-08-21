""" MultiQC module to parse output from ngs-disambiguate. """

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Disambiguate',
            anchor='disambiguate',
            href="https://github.com/AstraZeneca-NGS/disambiguate",
            info=("is an algorithm for disambiguating reads aligned "
                  "to two different species (e.g. human and mouse)."))

        self.data = {}

        for summary_file in self.find_log_files('disambiguate'):
            summ_data = self.parse_summary(summary_file['f'])

            for sample, counts in summ_data.items():
                # Clean sample name.
                sample = self.clean_s_name(sample, summary_file['root'])

                # Check for duplicates.
                if sample in self.data:
                    log.debug('Duplicate sample name found! Overwriting: %s',
                              sample)

                # Add to data and add data source.
                self.data[sample] = counts
                self.add_data_source(summary_file, s_name=sample)

        self.data = self.ignore_samples(self.data)

        if len(self.data) == 0:
            log.debug("Could not find any data in %s", config.analysis_dir)
            raise UserWarning

        log.info("Found %d reports", len(self.data))

        self.add_stats_table()
        self.add_stats_plot()
        self.write_stats_to_file()


    def parse_summary(self, contents):
        """Parses summary file into a dictionary of counts."""

        lines = contents.strip().split('\n')

        data = {}
        for row in lines[1:]:
            split = row.strip().split('\t')

            sample = split[0]

            data[sample] = {
                'species_a': int(split[1]),
                'species_b': int(split[2]),
                'ambiguous': int(split[3])
            }

        return data


    def add_stats_table(self):
        """Adds stats to general table."""

        totals = {sample: sum(counts.values())
                  for sample, counts in self.data.items()}

        percentages = {sample: {k: (v / totals[sample]) * 100
                                for k, v in counts.items()}
                       for sample, counts in self.data.items()}

        headers = {
            'species_a': {
                'title': '% Species a',
                'description': 'Percentage of reads mapping to species a',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'YlGn'
            }
        }

        self.general_stats_addcols(percentages, headers)


    def write_stats_to_file(self):
        """Writes stats to data file."""
        self.write_data_file(self.data, 'multiqc_disambiguate')


    def add_stats_plot(self):
        """Plots alignment stats as bargraph."""

        keys = OrderedDict()
        keys['species_a'] = {'color': '#437bb1', 'name': 'Species a'}
        keys['species_b'] = {'color': '#b1084c', 'name': 'Species b'}
        keys['ambiguous'] = {'color': '#333333', 'name': 'Ambiguous'}

        plot_config = {
            'id': "disambiguated_alignments",
            'title': "Disambiguated alignments",
            'cpswitch_counts_label': "# Reads",
            'ylab': "# Reads"
        }

        self.add_section(
            plot=bargraph.plot(self.data, keys, plot_config)
        )
