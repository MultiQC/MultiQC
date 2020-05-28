from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, beeswarm, table
from multiqc.modules.base_module import BaseMultiqcModule

import logging
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Malt Module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='MALT',
            anchor='malt',
            href='https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/malt/',
            info='performs alignment of metagenomic reads against a database of reference sequences (such as NR, GenBank or Silva) and produces a MEGAN RMA file as output.'
        )

        # Find and load Malt reports
        self.malt_raw_data = dict()

        for f in self.find_log_files('malt', filehandles=True):
            self.parse_logs(f)

        if len(self.malt_raw_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.malt_raw_data)))

        self.malt_summary_table()
        self.mappability_barplot()
        self.taxonomic_assignation_barplot()

    def parse_logs(self, f):
        """Parses a Malt log file"""
        data = {}
        mappability = {}
        taxo_success = {}
        reading = False
        for line in f['f']:
            line = line.rstrip()
            if line.startswith("+++++ Aligning file:") and reading == False:
                reading = True
                sample = line.split()[-1]
                sample = self.clean_s_name(sample, f['root'])
                data[sample] = {}
                continue
            if line.startswith("Total reads:") and reading:
                data[sample]["Total reads"] = int(
                    line.split()[-1].replace(",", ""))
                continue
            if line.startswith("Assig. Taxonomy:") and reading:
                data[sample]['Assig. Taxonomy'] = int(
                    line.split()[-1].replace(",", ""))
                continue
            if line.startswith("Num. of queries:") and reading:
                data[sample]["Num. of queries"] = int(
                    line.split()[-1].replace(",", ""))
                continue
            if line.startswith("Aligned queries:") and reading:
                data[sample]["Aligned queries"] = int(
                    line.split()[-1].replace(",", ""))
                continue
            if line.startswith("Num. alignments:") and reading:
                reading = False
                data[sample]['num. alignments'] = int(
                    line.split()[-1].replace(",", ""))
                continue
            self.malt_raw_data = data
            # self.general_stats_addcols(data)

    def mappability_barplot(self):
        """Mappability barplot

        Mappability = (Total reads / Num. of queries) * 100
        """
        data = self.malt_raw_data
        for samp in data:
            data[samp]['Non mapped'] = data[samp]["Num. of queries"] - \
                data[samp]['Total reads']
        cats = OrderedDict()
        cats['Total reads'] = {'name': 'Mapped reads'}
        cats['Non mapped'] = {'name': 'Non Mapped reads'}
        config = {
            'id': 'malt-mappability',
            'title': 'MALT: Metagenomic Mappability'
        }
        self.add_section(
            name='Metagenomic Mappability',
            anchor='malt-mappability',
            description='Number of mapped reads',
            helptext="""
                Mappability = (Assig. Taxonomy / Num. of queries) * 100
            """,
            plot=bargraph.plot(data, cats, config)
        )

    def taxonomic_assignation_barplot(self):
        """Taxonomic assignment barplot

        Taxonomic assignment success = (Assig. Taxonomy / Total reads) * 100
        """
        data = self.malt_raw_data
        for samp in data:
            data[samp]['No Assig. Taxonomy'] = data[samp]["Total reads"] - \
                data[samp]['Assig. Taxonomy']
        cats = ['Assig. Taxonomy', "No Assig. Taxonomy"]
        config = {
            'id': 'malt-taxonomic-success',
            'title': 'MALT: Taxonomic assignment success'
        }
        self.add_section(
            name='Taxonomic assignment success',
            anchor='malt-taxonomic-success',
            description='Number of mapped reads assigned to a taxonomic node',
            helptext="""
                Taxonomic assignment success = (Assig. Taxonomy / Total Reads) * 100
            """,
            plot=bargraph.plot(data, cats, config)
        )

    def malt_summary_table(self):
        """MALT summary table"""
        data = self.malt_raw_data
        for samp in data:
            data[samp]['Mappability'] = (
                data[samp]['Total reads'] / data[samp]["Num. of queries"])*100
            data[samp]['Taxonomic assignment success'] = (
                data[samp]['Assig. Taxonomy']/data[samp]['Total reads'])*100
        headers = OrderedDict()
        headers['Num. of queries'] = {
            'title': 'Number of reads in sample',
            'description': 'Number of reads in sample',
            'format': '{:d}'
        }
        headers['Total reads'] = {
            'title': 'Mapped reads',
            'description': 'Number of mapped reads',
            'format': '{:d}'
        }
        headers['Mappability'] = {
            'title': "% Metagenomic Mappability",
            'description': 'Percentage of mapped reads',
            'suffix': '%',
            'max': 100,
            'format': '{:,.2f}'
        }
        headers['Taxonomic assignment success'] = {
            'title': "% Taxonomic assignment success",
            'description': 'Percentage of mapped reads assigned to a taxonomic node',
            'suffix': '%',
            'max': 100,
            'format': '{:,.2f}'
        }
        headers['Assig. Taxonomy'] = {
            'title': 'Assig. Taxonomy',
            'description': 'Number of reads assigned to a Taxonomic node',
            'format': '{:d}',
            'hidden': True
        }
        config = {
            'namespace': 'malt'
        }
        self.add_section(
            name='MALT summary table',
            anchor='malt-summary-table',
            description='MALT summary statistics table',
            plot=table.plot(data, headers, config)
        )
        self.general_stats_addcols(data, headers)
