#!/usr/bin/env python

""" MultiQC module to parse output from Peddy """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Peddy module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Peddy', anchor='peddy',
        href='https://github.com/brentp/peddy',
        info="calculates genotype :: pedigree correspondence checks, ancestry checks and sex checks using VCF files.")

        # Find and load any Peddy reports
        self.peddy_data = dict()
        self.peddy_length_counts = dict()
        self.peddy_length_exp = dict()
        self.peddy_length_obsexp = dict()

        for f in self.find_log_files('peddy/summary_table'):
            parsed_data = self.parse_peddy_summary(f)
            if parsed_data is not None:
                for s_name in parsed_data:
                    s_name = self.clean_s_name(s_name, f['root'])
                    try:
                        self.peddy_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.peddy_data[s_name] = parsed_data[s_name]

        for pattern in ['het_check', 'ped_check', 'sex_check']:
            sp_key = 'peddy/{}'.format(pattern)
            for f in self.find_log_files(sp_key):
                # some columns have the same name in het_check and sex_check (median_depth)
                # pass pattern to parse_peddy_csv so the column names can include pattern to avoid being overwritten
                parsed_data = self.parse_peddy_csv(f, pattern)
                if parsed_data is not None:
                    for s_name in parsed_data:
                        s_name = self.clean_s_name(s_name, f['root'])
                        try:
                            self.peddy_data[s_name].update(parsed_data[s_name])
                        except KeyError:
                            self.peddy_data[s_name] = parsed_data[s_name]

        # Filter to strip out ignored sample names
        self.peddy_data = self.ignore_samples(self.peddy_data)

        if len(self.peddy_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.peddy_data)))

        # Write parsed report data to a file
        self.write_data_file(self.peddy_data, 'multiqc_peddy')

        # Basic Stats Table
        self.peddy_general_stats_table()

        # PCA plot
        self.peddy_pca_plot()

        # Relatedness plot
        self.peddy_relatedness_plot()

        # hetcheck plot
        self.peddy_het_check_plot()

    def parse_peddy_summary(self, f):
        """ Go through log file looking for peddy output """
        parsed_data = dict()
        headers = None
        for l in f['f'].splitlines():
            s = l.split("\t")
            if headers is None:
                s[0] = s[0].lstrip('#')
                headers = s
            else:
                parsed_data[s[1]] = dict()
                for i, v in enumerate(s):
                    if i != 1:
                        try:
                            parsed_data[s[1]][headers[i]] = float(v)
                        except ValueError:
                            parsed_data[s[1]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def parse_peddy_csv(self, f, pattern):
        """ Parse csv output from peddy """
        parsed_data = dict()
        headers = None
        s_name_idx = None
        for l in f['f'].splitlines():
            s = l.split(",")
            if headers is None:
                headers = s
                try:
                    s_name_idx = [headers.index("sample_id")]
                except ValueError:
                    try:
                        s_name_idx = [headers.index("sample_a"), headers.index("sample_b")]
                    except ValueError:
                        log.warn("Could not find sample name in Peddy output: {}".format(f['fn']))
                        return None
            else:
                s_name = '-'.join([s[idx] for idx in s_name_idx])
                parsed_data[s_name] = dict()
                for i, v in enumerate(s):
                    if i not in s_name_idx:
                        if headers[i] == "error" and pattern == "sex_check":
                            v = "True" if v == "False" else "False"
                        try:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = float(v)
                        except ValueError:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def peddy_general_stats_table(self):
        """ Take the parsed stats from the Peddy report and add it to the
        basic stats table at the top of the report """

        family_ids = [ x.get('family_id') for x in self.peddy_data.values() ]

        headers = OrderedDict()
        headers['family_id'] = {
            'title': 'Family ID',
            'hidden': True if all([v == family_ids[0] for v in family_ids]) else False
        }
        headers['ancestry-prediction'] = {
            'title': 'Ancestry',
            'description': 'Ancestry Prediction',
        }
        headers['sex_het_ratio'] = {
            'title': 'Sex / Het Ratio',
        }
        headers['error_sex_check'] = {
            'title': 'Correct Sex',
            'description': 'Displays False if error in sample sex prediction',
        }
        self.general_stats_addcols(self.peddy_data, headers)

    def peddy_pca_plot(self):

        data = dict()
        for s_name, d in self.peddy_data.items():
            if 'PC1_het_check' in d and 'PC2_het_check' in d:
                data[s_name] = {
                    'x': d['PC1_het_check'],
                    'y': d['PC2_het_check'],
                }

        pconfig = {
            'id': 'peddy_pca_plot',
            'title': 'Peddy: PCA Plot',
            'xlab': 'PC1',
            'ylab': 'PC2'
        }

        if len(data) > 0:
            self.add_section (
                name = 'PCA Plot',
                anchor = 'peddy-pca-plot',
                plot = scatter.plot(data, pconfig)
            )

    def peddy_relatedness_plot(self):
        data = dict()
        for s_name, d in self.peddy_data.items():
            if 'ibs0_ped_check' in d and 'ibs2_ped_check' in d:
                data[s_name] = {
                    'x': d['ibs0_ped_check'],
                    'y': d['ibs2_ped_check']
                }
            if 'rel_ped_check' in d:
                if d['rel_ped_check'] < 0.25:
                    data[s_name]['color'] = 'rgba(109, 164, 202, 0.9)'
                elif d['rel_ped_check'] < 0.5:
                    data[s_name]['color'] = 'rgba(250, 160, 81, 0.8)'
                else:
                    data[s_name]['color'] = 'rgba(43, 159, 43, 0.8)'

        pconfig = {
            'id': 'peddy_relatedness_plot',
            'title': 'Peddy: Relatedness Plot',
            'xlab': 'IBS0 (no alleles shared)',
            'ylab': 'IBS2 (both alleles shared)',
        }

        if len(data) > 0:
            self.add_section (
                name = 'Relatedness',
                anchor = 'peddy-relatedness-plot',
                description = """Shared allele rates between sample pairs. Points are coloured by degree of relatedness:
                <span style="color: #6DA4CA;">less than 0.25</span>,
                <span style="color: #FAA051;">0.25 - 0.5</span>,
                <span style="color: #2B9F2B;">greather than 0.5</span>.""",
                plot = scatter.plot(data, pconfig)
            )

    def peddy_het_check_plot(self):
        """plot the het_check scatter plot"""
        # empty dictionary to add sample names, and dictionary of values
        data = {}

        # for each sample, and list in self.peddy_data
        for s_name, d in self.peddy_data.items():
            # check the sample contains the required columns
            if 'median_depth_het_check' in d and 'het_ratio_het_check' in d:
                # add sample to dictionary with value as a dictionary of points to plot
                data[s_name] = {
                    'x': d['median_depth_het_check'],
                    'y': d['het_ratio_het_check']
                }

        pconfig = {
            'id': 'peddy_het_check_plot',
            'title': 'Peddy: Het Check',
            'xlab': 'median depth',
            'ylab': 'proportion het calls',
        }

        self.add_section (
            name = 'Het Check',
            description = "Proportion of sites that were heterozygous against median depth.",
            helptext = """
            A high proportion of heterozygous sites suggests contamination, a low proportion suggests consanguinity.
            
            See [the main peddy documentation](https://peddy.readthedocs.io/en/latest/output.html#het-check) for more details about the `het_check` command.
            """,
            anchor = 'peddy-hetcheck-plot',
            plot = scatter.plot(data, pconfig)
        )
