""" MultiQC module to parse output from fgbio """

from multiqc.modules.base_module import BaseMultiqcModule

from .groupreadsbyumi import GroupReadsByUmiMixin

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule, GroupReadsByUmiMixin):
    """ fgbio has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='fgbio',
            anchor='fgbio', target='fgbio',
            href='http://fulcrumgenomics.github.io/fgbio/',
            info=("  is a command line toolkit for working with genomic and "
                  "particularly next generation sequencing data.."))

        n = dict()
        n["groupreadsbyumi"] = self.parse_groupreadsbyumi_log()

        config = {
            'id': 'umi_support',
            'title': 'Familizy size count',
            'ylab': '# UMI',
            'xlab': 'Reads supporting UMI',
            "xmax": 15,
            'xDecimals': False
        }

        self.add_section(name = 'GroupReadsByUmi statistics',
        anchor = 'fgbio-groupreadsbyumi',
        description = 'During GroupReadsByUmi processing family size count data is generated, showing number of UMIs represented by a certain number of reads. ',
        helptext = '''**Note!** Multiqc expects the input file to have the following format <SAMPLE>.histo.tsv, SAMPLE will be for naming data.''',
        plot = linegraph.plot(n["groupreadsbyumi"], config))
