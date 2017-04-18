#!/usr/bin/env python

""" MultiQC module to parse output from Adapter Removal """

from __future__ import print_function
from collections import OrderedDict
import logging
import copy

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Adapter Removal',
        anchor='adapterRemoval', target='Adapter Removal',
        href='https://github.com/MikkelSchubert/adapterremoval',
        info=" rapid adapter trimming, identification, and read merging ")

        self.__read_type = None
        self.__collapsed = None
        self.s_name = None
        self.adapter_removal_data = {
            'single': dict(),
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }

        # variable definition for single- and paired-end reads
        self.arc_mate1 = {
            'single': dict(),
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }
        self.arc_discarged = {
            'single': dict(),
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }
        self.arc_all = {
            'single': dict(),
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }

        # variable definition for paired-end only reads
        self.arc_mate2 = {
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }
        self.arc_singleton = {
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }
        self.arc_collapsed = {
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }
        self.arc_collapsed_truncated = {
            'paired': {
                'collapsed': dict(),
                'noncollapsed': dict(),
            }
        }

        for f in self.find_log_files(config.sp['adapterRemoval'], filehandles=True):
            self.s_name = f['s_name']
            try:
                parsed_data = self.parse_settings_file(f)
            except UserWarning:
                continue
            if parsed_data is not None:
                if self.__read_type == 'single':
                    self.adapter_removal_data[self.__read_type][self.s_name] = parsed_data
                else:
                    if self.__collapsed:
                        self.adapter_removal_data[self.__read_type]['collapsed'][self.s_name] = parsed_data
                    else:
                        self.adapter_removal_data[self.__read_type]['noncollapsed'][self.s_name] = parsed_data

        if len(self.adapter_removal_data['single']) == 0 and len(self.adapter_removal_data['paired']) == 0:
            log.warning("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        elements_count = len(self.adapter_removal_data['single']) + \
                         len(self.adapter_removal_data['paired']['noncollapsed']) + \
                         len(self.adapter_removal_data['paired']['collapsed'])

        log.info("Found {} reports".format(elements_count))

        # Write parsed report data to a file
        self.write_data_file(self.adapter_removal_data['single'], 'multiqc_adapter_removal_single_end')
        self.write_data_file(self.adapter_removal_data['paired']['noncollapsed'], 'multiqc_adapter_removal_paired_end_noncollapsed')
        self.write_data_file(self.adapter_removal_data['paired']['collapsed'], 'multiqc_adapter_removal_paired_end_collapsed')

        # Start the sections
        self.sections = list()

        self.adapter_removal_counts_chart()
        self.adapter_removal_retained_chart()
        self.adapter_removal_length_dist_plot()

    def parse_settings_file(self, f):

        self.result_data = {
            'total': None,
            'unaligned': None,
            'aligned': None,
            'reads_total': None,
            'retained': None,
            'discarded': None,
        }

        settings_data = {'header': []}

        block_title = None
        for i, line in enumerate(f['f']):

            line = line.rstrip('\n')
            if line == '':
                continue

            if not block_title:
                block_title = 'header'
                settings_data[block_title].append(str(line))
                continue

            if line.startswith('['):
                block_title = str(line.strip('[]'))
                settings_data[block_title] = []
                continue

            settings_data[block_title].append(str(line))

        # set data for further working
        self.set_result_data(settings_data)

        return self.result_data

    def set_result_data(self, settings_data):
        # set read and collapsed type
        self.set_ar_type(settings_data['Length distribution'])

        # set result_data
        self.set_trim_stat(settings_data['Trimming statistics'])
        self.set_len_dist(settings_data['Length distribution'])

    def set_ar_type(self, len_dist_data):
        head_line = len_dist_data[0].rstrip('\n').split('\t')
        self.__read_type = 'paired' if head_line[2] == 'Mate2' else 'single'
        self.__collapsed = True if head_line[-3] == 'CollapsedTruncated' else False

        # biological/technical relevance is not clear -> skip
        if self.__read_type == 'single' and self.__collapsed:
            log.warning("Case single-end and collapse is not " \
                        "implemented -> File %s skipped" % self.s_name)
            raise UserWarning

    def set_trim_stat(self, trim_data):
        data_pattern = {'total': 0,
                        'unaligned': 1,
                        'aligned': 2,
                        'discarded_m1': 3,
                        'singleton_m1': 4,
                        'retained': 6}

        if self.__read_type == 'paired':
            data_pattern['discarded_m2'] = 5
            data_pattern['singleton_m2'] = 6
            if not self.__collapsed:
                data_pattern['retained'] = 8
            else:
                data_pattern['full-length_cp'] = 8
                data_pattern['truncated_cp'] = 9
                data_pattern['retained'] = 10

        for title, key in data_pattern.items():
            tmp = trim_data[key]
            value = tmp.split(': ')[1]
            self.result_data[title] = int(value)

        reads_total = self.result_data['total']
        aligned_total = self.result_data['aligned']
        unaligned_total = self.result_data['unaligned']
        if self.__read_type == 'paired':
            reads_total = self.result_data['total'] * 2
            aligned_total = self.result_data['aligned'] * 2
            unaligned_total = self.result_data['unaligned'] * 2

        self.result_data['aligned_total'] = aligned_total
        self.result_data['unaligned_total'] = unaligned_total
        self.result_data['reads_total'] = reads_total
        self.result_data['discarded_total'] = reads_total - self.result_data['retained']

        if self.__read_type == 'paired':
            self.result_data['retained_reads'] = self.result_data['retained'] - self.result_data['singleton_m1'] - self.result_data['singleton_m2']

    def set_len_dist(self, len_dist_data):

        for line in len_dist_data[1:]:
            l_data = line.rstrip('\n').split('\t')
            l_data = map(int, l_data)
            if self.__read_type == 'single':
                if not self.__collapsed:
                    if self.s_name not in self.arc_mate1['single']:
                        self.arc_mate1['single'][self.s_name] = dict()
                    self.arc_mate1['single'][self.s_name][l_data[0]] = l_data[1]

                    if self.s_name not in self.arc_discarged['single']:
                        self.arc_discarged['single'][self.s_name] = dict()
                    self.arc_discarged['single'][self.s_name][l_data[0]] = l_data[2]

                    if self.s_name not in self.arc_all['single']:
                        self.arc_all['single'][self.s_name] = dict()
                    self.arc_all['single'][self.s_name][l_data[0]] = l_data[3]
                else:
                    # this case should not be reached (see case at method set_ar_type())
                    pass
            else:
                if not self.__collapsed:
                    if self.s_name not in self.arc_mate1['paired']['noncollapsed']:
                        self.arc_mate1['paired']['noncollapsed'][self.s_name] = dict()
                    self.arc_mate1['paired']['noncollapsed'][self.s_name][l_data[0]] = l_data[1]

                    if self.s_name not in self.arc_mate2['paired']['noncollapsed']:
                        self.arc_mate2['paired']['noncollapsed'][self.s_name] = dict()
                    self.arc_mate2['paired']['noncollapsed'][self.s_name][l_data[0]] = l_data[2]

                    if self.s_name not in self.arc_singleton['paired']['noncollapsed']:
                        self.arc_singleton['paired']['noncollapsed'][self.s_name] = dict()
                    self.arc_singleton['paired']['noncollapsed'][self.s_name][l_data[0]] = l_data[3]

                    if self.s_name not in self.arc_discarged['paired']['noncollapsed']:
                        self.arc_discarged['paired']['noncollapsed'][self.s_name] = dict()
                    self.arc_discarged['paired']['noncollapsed'][self.s_name][l_data[0]] = l_data[4]

                    if self.s_name not in self.arc_all['paired']['noncollapsed']:
                        self.arc_all['paired']['noncollapsed'][self.s_name] = dict()
                    self.arc_all['paired']['noncollapsed'][self.s_name][l_data[0]] = l_data[5]
                else:
                    if self.s_name not in self.arc_mate1['paired']['collapsed']:
                        self.arc_mate1['paired']['collapsed'][self.s_name] = dict()
                    self.arc_mate1['paired']['collapsed'][self.s_name][l_data[0]] = l_data[1]

                    if self.s_name not in self.arc_mate2['paired']['collapsed']:
                        self.arc_mate2['paired']['collapsed'][self.s_name] = dict()
                    self.arc_mate2['paired']['collapsed'][self.s_name][l_data[0]] = l_data[2]

                    if self.s_name not in self.arc_singleton['paired']['collapsed']:
                        self.arc_singleton['paired']['collapsed'][self.s_name] = dict()
                    self.arc_singleton['paired']['collapsed'][self.s_name][l_data[0]] = l_data[3]

                    if self.s_name not in self.arc_collapsed['paired']['collapsed']:
                        self.arc_collapsed['paired']['collapsed'][self.s_name] = dict()
                    self.arc_collapsed['paired']['collapsed'][self.s_name][l_data[0]] = l_data[4]

                    if self.s_name not in self.arc_collapsed_truncated['paired']['collapsed']:
                        self.arc_collapsed_truncated['paired']['collapsed'][self.s_name] = dict()
                    self.arc_collapsed_truncated['paired']['collapsed'][self.s_name][l_data[0]] = l_data[5]

                    if self.s_name not in self.arc_discarged['paired']['collapsed']:
                        self.arc_discarged['paired']['collapsed'][self.s_name] = dict()
                    self.arc_discarged['paired']['collapsed'][self.s_name][l_data[0]] = l_data[6]

                    if self.s_name not in self.arc_all['paired']['collapsed']:
                        self.arc_all['paired']['collapsed'][self.s_name] = dict()
                    self.arc_all['paired']['collapsed'][self.s_name][l_data[0]] = l_data[7]

    def adapter_removal_counts_chart(self):

        cats = OrderedDict()
        cats['aligned_total'] = {'name': 'with adapter'}
        cats['unaligned_total'] = {'name': 'without adapter'}
        pconfig = {
            'title': 'Adapter Alignments',
            'ylab': '# Reads',
            'hide_zero_cats': False,
            'cpswitch_counts_label': 'Number of Reads'
        }

        # plot different results if exists
        if self.adapter_removal_data['single']:
            pconfig['id'] = 'ar_alignment_plot_se'
            self.sections.append({
                'name': 'Adapter Alignments Single-End',
                'anchor': 'ar_alignment_se',
                'content': '<p>The proportions of reads with and without adapter.</p>' +
                           bargraph.plot(self.adapter_removal_data['single'], cats, pconfig)
            })

        if self.adapter_removal_data['paired']['noncollapsed']:
            pconfig['id'] = 'ar_alignment_plot_penc'
            self.sections.append({
                'name': 'Adapter Alignments Paired-End Noncollapsed',
                'anchor': 'adapter_removal_alignment_penc',
                'content': '<p>The proportions of reads with and without adapter.</p>' +
                           bargraph.plot(self.adapter_removal_data['paired']['noncollapsed'], cats, pconfig)
            })

        if self.adapter_removal_data['paired']['collapsed']:
            pconfig['id'] = 'ar_alignment_plot_pec'
            self.sections.append({
                'name': 'Adapter Alignments Paired-End Collapsed',
                'anchor': 'adapter_removal_alignment_pec',
                'content': '<p>The proportions of reads with and without adapter.</p>' +
                           bargraph.plot(self.adapter_removal_data['paired']['collapsed'], cats, pconfig)
            })

    def adapter_removal_retained_chart(self):

        pconfig = {
            'title': 'retained and discarded',
            'ylab': '# Reads',
            'hide_zero_cats': False,
            'cpswitch_counts_label': 'Number of Reads'
        }

        # plot different results if exists
        if self.adapter_removal_data['single']:
            cats_se = OrderedDict()
            cats_se['singleton_m1'] = {
                'name': 'singleton mate1',
                'color': '#305ead'
            }
            cats_se['discarded_m1'] = {
                'name': 'discarded mate1',
                'color': '#bd5dbf'
            }
            pconfig['id'] = 'ar_retained_plot_se'
            self.sections.append({
                'name': 'Retained and Discarded Single-End',
                'anchor': 'adapter_removal_retained_plot_se',
                'content': '<p>The proportions of retained and discarded reads.</p>' +
                           bargraph.plot(self.adapter_removal_data['single'], cats_se, pconfig)
            })

        if self.adapter_removal_data['paired']['noncollapsed']:
            cats_penc = OrderedDict()
            cats_penc['singleton_m1'] = {
                'name': 'singleton mate1',
                'color': '#305ead'
            }
            cats_penc['singleton_m2'] = {
                'name': 'singleton mate2',
                'color': '#3b66af'
            }
            cats_penc['retained_reads'] = {
                'name': 'retained read pairs',
                'color': '#b1b5bc'
            }
            cats_penc['discarded_m1'] = {
                'name': 'discarded mate1',
                'color': '#e08fc1'
            }
            cats_penc['discarded_m2'] = {
                'name': 'discarded mate2',
                'color': '#bd5dbf'
            }
            pconfig['id'] = 'ar_retained_plot_penc'
            self.sections.append({
                'name': 'Retained and Discarded Paired-End Noncollapsed',
                'anchor': 'adapter_removal_retained_plot_penc',
                'content': '<p>The proportions of retained and discarded reads.</p>' +
                           bargraph.plot(self.adapter_removal_data['paired']['noncollapsed'], cats_penc, pconfig)
            })

        if self.adapter_removal_data['paired']['collapsed']:
            cats_pec = OrderedDict()
            cats_pec['singleton_m1'] = {
                'name': 'singleton mate1',
                'color': '#305ead'
            }
            cats_pec['singleton_m2'] = {
                'name': 'singleton mate2',
                'color': '#3b66af'
            }
            cats_pec['retained_reads'] = {
                'name': 'retained read pairs',
                'color': '#b1b5bc'
            }
            cats_pec['full-length_cp'] = {
                'name': 'full-length collapsed pairs',
                'color': '#717f99'
            }
            cats_pec['truncated_cp'] = {
                'name': 'truncated collapsed pairs',
                'color': '#748ab2'
            }
            cats_pec['discarded_m1'] = {
                'name': 'discarded mate1',
                'color': '#e08fc1'
            }
            cats_pec['discarded_m2'] = {
                'name': 'discarded mate2',
                'color': '#bd5dbf'
            }
            pconfig['id'] = 'ar_retained_plot_pec'
            self.sections.append({
                'name': 'Retained and Discarded Paired-End Collapsed',
                'anchor': 'adapter_removal_retained_plot_pec',
                'content': '<p>The proportions of retained and discarded reads.</p>' +
                           bargraph.plot(self.adapter_removal_data['paired']['collapsed'], cats_pec, pconfig)
            })

    def adapter_removal_length_dist_plot(self):

        config_template = {
            'title': 'Length Distribution',
            'ylab': 'Counts',
            'xlab': 'read length',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': []
        }

        pconfig = {
            'single': copy.deepcopy(config_template),
            'paired': {
                'collapsed': copy.deepcopy(config_template),
                'noncollapsed': copy.deepcopy(config_template),
            }
        }

        dl_mate1 = {'name': 'Mate1', 'ylab': 'Count'}
        dl_mate2 = {'name': 'Mate2', 'ylab': 'Count'}
        dl_singleton = {'name': 'Singleton', 'ylab': 'Count'}
        dl_collapsed = {'name': 'Collapsed', 'ylab': 'Count'}
        dl_collapsed_truncated = {'name': 'Collapsed Truncated', 'ylab': 'Count'}
        dl_discarded = {'name': 'Discarded', 'ylab': 'Count'}
        dl_all = {'name': 'All', 'ylab': 'Count'}

        lineplot_data = {
            'single': None,
            'paired': {
                'collapsed': None,
                'noncollapsed': None,
            }
        }

        if self.adapter_removal_data['single']:
            lineplot_data['single'] = [
                self.arc_mate1['single'],
                self.arc_discarged['single'],
                self.arc_all['single']]
            pconfig['single']['id'] = 'ar_lenght_count_plot_se'
            pconfig['single']['data_labels'].extend([
                dl_mate1,
                dl_discarded,
                dl_all])
            self.sections.append({
                'name': 'Lenght Distribution Single End',
                'anchor': 'ar_lenght_count_se',
                'content': '<p>The lenght distribution of reads after processing adapter alignment.</p>' +
                           linegraph.plot(lineplot_data['single'], pconfig['single'])
            })

        if self.adapter_removal_data['paired']['noncollapsed']:
            lineplot_data['paired']['noncollapsed'] = [
                self.arc_mate1['paired']['noncollapsed'],
                self.arc_mate2['paired']['noncollapsed'],
                self.arc_singleton['paired']['noncollapsed'],
                self.arc_discarged['paired']['noncollapsed'],
                self.arc_all['paired']['noncollapsed']]
            pconfig['paired']['noncollapsed']['id'] = 'ar_lenght_count_plot_penc'
            pconfig['paired']['noncollapsed']['data_labels'].extend([
                dl_mate1,
                dl_mate2,
                dl_singleton,
                dl_discarded,
                dl_all])
            self.sections.append({
                'name': 'Lenght Distribution Paired End Noncollapsed',
                'anchor': 'ar_lenght_count_penc',
                'content': '<p>The lenght distribution of reads after processing adapter alignment.</p>' +
                           linegraph.plot(lineplot_data['paired']['noncollapsed'], pconfig['paired']['noncollapsed'])
            })

        if self.adapter_removal_data['paired']['collapsed']:
            lineplot_data['paired']['collapsed'] = [
                self.arc_mate1['paired']['collapsed'],
                self.arc_mate2['paired']['collapsed'],
                self.arc_singleton['paired']['collapsed'],
                self.arc_collapsed['paired']['collapsed'],
                self.arc_collapsed_truncated['paired']['collapsed'],
                self.arc_discarged['paired']['collapsed'],
                self.arc_all['paired']['collapsed']]
            pconfig['paired']['collapsed']['id'] = 'ar_lenght_count_plot_pec'
            pconfig['paired']['collapsed']['data_labels'].extend([
                dl_mate1,
                dl_mate2,
                dl_singleton,
                dl_collapsed,
                dl_collapsed_truncated,
                dl_discarded,
                dl_all])
            self.sections.append({
                'name': 'Lenght Distribution Paired End Collapsed',
                'anchor': 'ar_lenght_count_pec',
                'content': '<p>The lenght distribution of reads after processing adapter alignment.</p>' +
                           linegraph.plot(lineplot_data['paired']['collapsed'], pconfig['paired']['collapsed'])
            })