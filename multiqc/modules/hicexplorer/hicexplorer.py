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
            self.mod_data[file['s_name']]['File'][0] = self.clean_s_name(self.mod_data[file['s_name']]['File'][0], file['root'])
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

        # key lists for plotting
        keys_list_plot_pairs = ['Pairs used', 'One mate unmapped', 'One mate not unique', 'One mate low quality',
                                'Dangling end', 'Self ligation (removed)', 'One mate not close to rest site', 'Same fragment',
                                'Self circle', 'Duplicated pairs']
        keys_list_contact_distance = ['Inter chromosomal', 'Short range', 'Long range']
        keys_list_read_orientation = ['Inward pairs', 'Outward pairs', 'Left pairs', 'Right pairs']

        # prepare the detail report section
        self.add_section(
            name='Pair categorization',
            anchor='hicexplorer_pairs_categorized',
            plot=self.hicexplorer_create_plot(keys_list_plot_pairs, 'HiCExplorer: Categorization of reads'),
            description='This figure contains the number of reads that were finally used to build the '
            'Hi-C matrix along with the reads that where filtered out.',
            helptext='''        
            #### Dangling ends

            These are reads that start with the restriction site and constitute reads that were digested but no ligated.

            #### Same fragment

            These are read mates, facing inward, separated by up to 800 bp that do not have a restriction enzyme in between.
            These read pairs are not valid Hi-C pairs.

            #### Self circle

            Self circles are defined as pairs within 25kb with \'outward\' read orientation

            #### Self ligation

            These are read pairs with a restriction site in between that are within 800 bp.'''
        )

        self.add_section(
            name='Contact distance',
            anchor='hicexplorer_contact_distance',
            plot=self.hicexplorer_create_plot(keys_list_contact_distance, 'Contact distance'),
            description='This figure contains information about the distance and location of the valid pairs used.',
            helptext='''  
            ####Long range
            
            Pairs with a distance greater than 20 kilobases

            ####Short range
            
            Pairs with a distance less than 20 kilobases

            ####Inter chromosomal

            Interchromosomal pairs.
            '''
        )

        self.add_section(
            name='Read orientation',
            anchor='hicexplorer_read_orientation',
            plot=self.hicexplorer_create_plot(keys_list_read_orientation, 'Read orientation'),
            description='This figure contains information about the orientation of the read pairs.',
            helptext='''  
                ####Inward pairs
                
                First mate is a forward read, second is reverse.
                --------------->              <----------------
                
                ####Outward pairs
                
                First mate is a reverse read, second is forward.
                --------------->              <----------------
                
                ####Left pairs
                
                Both are reverse reads.
                <---------------              <----------------
                
                ####Right pairs
                
                Both are forward reads.
                --------------->              ---------------->

            '''
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
                except ValueError:
                    data_.append(i)
            if len(data_) == 0:
                continue
            if s[0].startswith('short range'):
                s[0] = 'short range'
            elif s[0].startswith('same fragment'):
                s[0] = 'same fragment'
            s[0] = s[0].capitalize()
            data[s[0]] = data_
        return data

    def hicexplorer_basic_statistics(self):
        """Create the general statistics for HiCExplorer."""
        data = {}

        for file in self.mod_data:
            total_pairs = self.mod_data[file]['Pairs considered'][0]
            data_ = {
                'Pairs considered': self.mod_data[file]['Pairs considered'][0],

                'Pairs used': self.mod_data[file]['Pairs used'][0] / total_pairs,
                'Mapped': self.mod_data[file]['One mate unmapped'][0] / total_pairs,
                'Min rest. site distance': self.mod_data[file]['Min rest. site distance'][0],
                'Max rest. site distance': self.mod_data[file]['Max rest. site distance'][0],
            }
            data[self.mod_data[file]['File'][0]] = data_
        headers = OrderedDict()

        headers['Pairs considered'] = {
            'title': '{} Pairs'.format(config.read_count_prefix),
            'shared_key': 'read_count'
        }

        headers['Pairs used'] = {
            'title': '% Used pairs',
            'max': 100,
            'min': 0,
            'modify': lambda x: x * 100,
            'suffix': '%'
        }
        headers['Mapped'] = {
            'title': '% Mapped',
            'max': 100,
            'min': 0,
            'modify': lambda x: (1 - x) * 100,
            'scale': 'RdYlGn',
            'suffix': '%'
        }

        headers['Min rest. site distance'] = {
            'title': 'Min RE dist',
            'descripton': 'Minimum restriction site distance (bp)',
            'scale': None,
            'modify': lambda x: int(x)
        }

        headers['Max rest. site distance'] = {
            'title': 'Max RE dist',
            'descripton': 'Maximum restriction site distance (bp)',
            'scale': None,
            'modify': lambda x: int(x)
        }
        self.general_stats_addcols(data, headers)

    def hicexplorer_create_plot(self, pKeyList, pTitle):
        """Create the graphics containing information about the read quality."""

        keys = OrderedDict()

        for i, key_ in enumerate(pKeyList):
            keys[key_] = {'color': self.colors[i]}

        data = {}

        for data_ in self.mod_data:
            data['{}'.format(self.mod_data[data_]['File'][0])] = {}
            for key_ in pKeyList:
                data['{}'.format(self.mod_data[data_]['File'][0])][key_] = self.mod_data[data_][key_][0]

        config = {
            'title': pTitle,
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(data, keys, config)
