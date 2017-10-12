from multiqc.modules.base_module import BaseMultiqcModule
import logging
from collections import OrderedDict

from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HiCExplorer', anchor='hicexplorer',
                                            href="https://hicexplorer.readthedocs.io",
                                            info=" addresses the common tasks of Hi-C analysis from processing to visualization.")

        self.mod_data = dict()
        for file in self.find_log_files('hicexplorer'):
            self.mod_data[file['s_name']] = self.parse_logs(file['f'])
            self.add_data_source(file)
        if len(self.mod_data) == 0:
            raise UserWarning
        self.colors = ["#1f77b4",
                       "#ff7f0e",
                       "#2ca02c",
                       "#d62728",
                       "#9467bd",
                       "#8c564b",
                       "#7f7f7f",
                       "#bcbd22",
                       "#17becf",
                       "#D2691E"]
        # prepare the basic statistics for hicexplorer
        self.hicexplorer_basic_statistics()

        # prepare the detail report section
        self.add_section(
            name='Pair categorization',
            anchor='hicexplorer_pairs_categorized',
            plot=self.hicexplorer_plot_pairs_categorized(),
            description='This figure contains the number of reads that were finally used to build the '
                        'Hi-C matrix along with the reads that where filtered out.'
                        '<br><br><b>Dangling ends</b><br>'
                        'These are reads that start with the restriction site and constitute reads that were digested but no ligated.'
                        '<br><br><b>Same fragment</b><br>'
                        'These are read mates, facing inward, separated by up to 800 bp that do not have a restriction enzyme in between.'
                        ' These read pairs are not valid Hi-C pairs.'
                        '<br><br><b>Self circle</b><br>'
                        'Self circles are defined as pairs within 25kb with \'outward\' read orientation'
                        '<br><br><b>Self ligation</b><br>'
                        'These are read pairs with a restriction site in between that are within 800 bp.'
        )

        self.add_section(
            name='Contact distance',
            anchor='hicexplorer_contact_distance',
            plot=self.hicexplorer_contact_distance()
        )

        self.add_section(
            name='Read orientation',
            anchor='hicexplorer_read_orientation',
            plot=self.hicexplorer_read_orientation()
        )

    def parse_logs(self, f):
        """Parse a given HiCExplorer log file from hicBuildMatrix."""
        data = {}
        for l in f.splitlines():
            # catch empty lines
            if len(l) == 0:
                continue
            s = l.split("\t")
            data_ = []
            # catch lines with descriptive content: "Of pairs used:"
            for i in s[1:]:
                if len(i) == 0:
                    continue
                try:
                    i.replace('(', '')
                    i.replace(')', '')
                    i.replace(',', '')

                    data_.append(float(i))
                except:
                    data_.append(i)
            if len(data_) == 0:
                continue
            if s[0].startswith('short range'):
                s[0] = 'short range'
            elif s[0].startswith('same fragment'):
                s[0] = 'same fragment'
            data[s[0]] = tuple(data_)
        return data

    def hicexplorer_basic_statistics(self):
        """Create the general statistics for HiCExplorer."""
        data = {}

        for file in self.mod_data:
            total_pairs = self.mod_data[file]['Pairs considered'][0]
            data_ = {
                'Pairs considered': self.mod_data[file]['Pairs considered'][0],

                'Pairs used': self.mod_data[file]['Pairs used'][0] / total_pairs,
                'Unmapped': self.mod_data[file]['One mate unmapped'][0] / total_pairs,
                'Min rest. site distance': self.mod_data[file]['Min rest. site distance'][0],
                'Max rest. site distance': self.mod_data[file]['Max rest. site distance'][0],
            }
            data[self.mod_data[file]['File'][0]] = data_
        headers = OrderedDict()

        headers['Pairs considered'] = {
            'title': 'Number of pairs',
            'modify': lambda x: int(x)
        }

        headers['Pairs used'] = {
            'title': '% Used pairs',
            'max': 100,
            'min': 0,
            'modify': lambda x: x * 100
        }
        headers['Unmapped'] = {
            'title': '% Unmapped',
            'max': 100,
            'min': 0,
            'modify': lambda x: x * 100
        }

        headers['Min rest. site distance'] = {
            'title': 'Min rest. site distance',
            'scale': None,
        }

        headers['Max rest. site distance'] = {
            'title': 'Max rest. site distance',
            'scale': None,
        }
        self.general_stats_addcols(data, headers)

    def hicexplorer_plot_pairs_categorized(self):
        """Create the graphics containing information about the categorization of the used pairs."""

        keys = OrderedDict()
        keys['Pairs used'] = {'color': self.colors[0], 'name': 'Pairs used'}
        keys['One mate unmapped'] = {'color': self.colors[1], 'name': 'One mate unmapped'}
        keys['One mate not unique'] = {'color': self.colors[2], 'name': 'One mate not unique'}
        keys['One mate low quality'] = {'color': self.colors[3], 'name': 'One mate low quality'}
        keys['Dangling end'] = {'color': self.colors[4], 'name': 'Dangling end'}

        keys['Self ligation'] = {'color': self.colors[5], 'name': 'Self ligation'}
        keys['One mate not close to rest site'] = {'color': self.colors[6], 'name': 'One mate not close to rest site'}
        keys['Same fragment'] = {'color': self.colors[7], 'name': 'Same fragment'}
        keys['Self circle'] = {'color': self.colors[8], 'name': 'Self circle'}
        keys['Duplicated pairs'] = {'color': self.colors[9], 'name': 'Duplicated pairs'}

        data = {}

        for data_ in self.mod_data:
            data['{}'.format(self.mod_data[data_]['File'][0])] = {}
            data['{}'.format(self.mod_data[data_]['File'][0])]['Pairs used'] = self.mod_data[data_]['Pairs used'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['One mate unmapped'] = self.mod_data[data_]['One mate unmapped'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['One mate not unique'] = self.mod_data[data_]['One mate not unique'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['One mate low quality'] = self.mod_data[data_]['One mate low quality'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['Dangling end'] = self.mod_data[data_]['dangling end'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['Self ligation'] = self.mod_data[data_]['self ligation (removed)'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['One mate not close to rest site'] = self.mod_data[data_]['One mate not close to rest site'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['Same fragment'] = self.mod_data[data_]['same fragment'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['Self circle'] = self.mod_data[data_]['self circle'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['Duplicated pairs'] = self.mod_data[data_]['duplicated pairs'][0]

        config = {
            'title': 'HiCExplorer: Categorization of reads',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)

    def hicexplorer_contact_distance(self):
        """Create the plot for the contact distance information of the aligned reads."""

        keys = OrderedDict()
        keys['inter chromosomal'] = {'color': self.colors[0], 'name': 'Inter chromosomal'}
        keys['short range < 20kb'] = {'color': self.colors[1], 'name': 'Short range < 20kb'}
        keys['long range'] = {'color': self.colors[2], 'name': 'Long range'}

        data = {}

        for data_ in self.mod_data:
            data['{}'.format(self.mod_data[data_]['File'][0])] = {}

            data['{}'.format(self.mod_data[data_]['File'][0])]['inter chromosomal'] = self.mod_data[data_]['inter chromosomal'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['short range < 20kb'] = self.mod_data[data_]['short range'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['long range'] = self.mod_data[data_]['long range'][0]

        config = {
            'title': 'HiCExplorer: Contact distance',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)

    def hicexplorer_read_orientation(self):
        """Create the plot about the read orientation of the aligned reads."""

        keys = OrderedDict()
        keys['inward pairs'] = {'color': self.colors[0], 'name': 'Inward pairs'}
        keys['outward pairs'] = {'color': self.colors[1], 'name': 'Outward pairs'}
        keys['left pairs'] = {'color': self.colors[2], 'name': 'Left pairs'}
        keys['right pairs'] = {'color': self.colors[3], 'name': 'Right pairs'}

        data = {}

        for data_ in self.mod_data:
            data['{}'.format(self.mod_data[data_]['File'][0])] = {}

            data['{}'.format(self.mod_data[data_]['File'][0])]['inward pairs'] = self.mod_data[data_]['inward pairs'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['outward pairs'] = self.mod_data[data_]['outward pairs'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['left pairs'] = self.mod_data[data_]['left pairs'][0]
            data['{}'.format(self.mod_data[data_]['File'][0])]['right pairs'] = self.mod_data[data_]['right pairs'][0]

        config = {
            'title': 'HiCExplorer: Read orientation',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)
