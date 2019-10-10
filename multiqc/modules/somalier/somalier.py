#!/usr/bin/env python

""" MultiQC module to parse output from somalier """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

import random
from multiqc import config
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    somalier module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='somalier', anchor='somalier',
        href='https://github.com/brentp/somalier',
        info="calculates genotype :: pedigree correspondence checks from sketches derived from BAM/CRAM or VCF")

        # Find and load any somalier reports
        self.somalier_data = dict()
        self.somalier_length_counts = dict()
        self.somalier_length_exp = dict()
        self.somalier_length_obsexp = dict()

        # parse somalier summary file
        for f in self.find_log_files('somalier/samples'):
            parsed_data = self.parse_somalier_summary(f)
            if parsed_data is not None:
                for s_name in parsed_data:
                    s_name = self.clean_s_name(s_name, f['root'])
                    try:
                        self.somalier_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.somalier_data[s_name] = parsed_data[s_name]
        # parse somalier CSV files
        for pattern in ['pairs']:
            sp_key = 'somalier/{}'.format(pattern)
            for f in self.find_log_files(sp_key):
                # some columns have the same name in het_check and sex_check (median_depth)
                # pass pattern to parse_somalier_csv so the column names can include pattern to
                # avoid being overwritten
                parsed_data = self.parse_somalier_csv(f, pattern)
                if parsed_data is not None:
                    for s_name in parsed_data:
                        s_name = self.clean_s_name(s_name, f['root'])
                        try:
                            self.somalier_data[s_name].update(parsed_data[s_name])
                        except KeyError:
                            self.somalier_data[s_name] = parsed_data[s_name]

        # Filter to strip out ignored sample names
        self.somalier_data = self.ignore_samples(self.somalier_data)

        if len(self.somalier_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.somalier_data)))

        # Write parsed report data to a file
        self.write_data_file(self.somalier_data, 'multiqc_somalier')

        # Basic Stats Table
        self.somalier_general_stats_table()

        # Relatedness plot
        self.somalier_relatedness_plot()

        # hetcheck plot
        self.somalier_het_check_plot()

        self.somalier_sex_check_plot()

    def parse_somalier_summary(self, f):
        """ Go through log file looking for somalier output """
        parsed_data = dict()
        headers = None
        sample_i = -100
        for l in f['f'].splitlines():
            s = l.split("\t")
            if headers is None:
                s[0] = s[0].lstrip('#')
                headers = s
                sample_i = headers.index("sample")
            else:
                parsed_data[s[sample_i]] = dict()
                for i, v in enumerate(s):
                    if i != sample_i:
                        try:
                            parsed_data[s[sample_i]][headers[i]] = float(v)
                        except ValueError:
                            parsed_data[s[sample_i]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def parse_somalier_csv(self, f, pattern):
        """ Parse csv output from somalier """
        parsed_data = dict()
        headers = None
        s_name_idx = None
        for l in f['f'].splitlines():
            s = l.lstrip('#').split(",")
            if headers is None:
                headers = s
                try:
                    s_name_idx = [headers.index("sample")]
                except ValueError:
                    try:
                        s_name_idx = [headers.index("sample_a"), headers.index("sample_b")]
                    except ValueError:
                        log.warn("Could not find sample name in somalier output: {}".format(f['fn']))
                        return None
            else:
                s_name = '-'.join([s[idx] for idx in s_name_idx])
                parsed_data[s_name] = dict()
                for i, v in enumerate(s):
                    if i not in s_name_idx:
                        try:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = float(v)
                        except ValueError:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def somalier_general_stats_table(self):
        """ Take the parsed stats from the somalier report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['X_het'] = {
            'title': 'heterozygous variants on X chromosome',
        }
        headers['X_hom_alt'] = {
            'title': 'homozygous alternate variants on X chromosome',
        }
        self.general_stats_addcols(self.somalier_data, headers)

    def somalier_relatedness_plot(self):
        data = dict()
        for s_name, d in self.somalier_data.items():
            if 'ibs0_pairs' in d and 'ibs2_pairs' in d:
                data[s_name] = {
                    'x': d['ibs0_pairs'],
                    'y': d['ibs2_pairs']
                }
            if 'relatedness_pairs' in d:
                if d['expected_relatedness_pairs'] < 0.125:
                    data[s_name]['color'] = 'rgba(109, 164, 202, 0.9)'
                elif d['expected_relatedness_pairs'] < 0.5:
                    data[s_name]['color'] = 'rgba(250, 160, 81, 0.8)'
                else:
                    data[s_name]['color'] = 'rgba(43, 159, 43, 0.8)'

        pconfig = {
            'id': 'somalier_relatedness_plot',
            'title': 'somalier: Relatedness Plot',
            'xlab': 'IBS0 (no alleles shared)',
            'ylab': 'IBS2 (both alleles shared)',
        }

        if len(data) > 0:
            self.add_section (
                name = 'Relatedness',
                anchor = 'somalier-relatedness-plot',
                description = """Shared allele rates between sample pairs.  Points are coloured by degree of expected-relatedness:
                <span style="color: #6DA4CA;">less than 0.25</span>,
                <span style="color: #FAA051;">0.25 - 0.5</span>,
                <span style="color: #2B9F2B;">greather than 0.5</span>.""",
                plot = scatter.plot(data, pconfig)
            )

    def somalier_het_check_plot(self):
        """plot the het_check scatter plot"""
        # empty dictionary to add sample names, and dictionary of values
        data = {}

        # for each sample, and list in self.somalier_data
        for s_name, d in self.somalier_data.items():
            # check the sample contains the required columns
            if 'gt_depth_mean' in d and 'ab_std' in d:
                # add sample to dictionary with value as a dictionary of points to plot
                data[s_name] = {
                    'x': d['gt_depth_mean'],
                    'y': d['ab_std']
                }

        pconfig = {
            'id': 'somalier_het_check_plot',
            'title': 'somalier: Het Check',
            'xlab': 'mean depth',
            'ylab': 'standard deviation of allele-balance',
        }

        self.add_section (
            name = 'Het Check',
            description = "Std devation of heterozygous allele balance against mean depth.",
            helptext = """A high standard deviation in allele balance suggests contamination.
            """,
            anchor = 'somalier-hetcheck-plot',
            plot = scatter.plot(data, pconfig)
        )

    def somalier_sex_check_plot(self):
        data = {}
        sex_index = {"female": 0, "male": 1, "unknown": 2}

        for s_name, d in self.somalier_data.items():
            if 'X_depth_mean' in d and 'pedigree_sex' in d:
                data[s_name] = {
                    'x': (random.random() - 0.5) * 0.1 + sex_index.get(d['pedigree_sex'], 2),
                    'y': d["X_depth_mean"]
                }

        pconfig = {
            'id': 'somalier_sex_check_plot',
            'title': 'somalier: Sex Check',
            'xlab': 'Sex From Ped',
            'ylab': 'Scaled mean depth on X',
            'categories': ["Female", "Male", "Unknown"]
        }

        self.add_section(
            name = 'Sex Check',
            description = "Predicted sex against scaled depth on X",
            helptext = """
            Higher values of depth, low values suggest male.
            """,
            anchor='somalier-sexcheck-plot',
            plot=scatter.plot(data, pconfig)
        )
