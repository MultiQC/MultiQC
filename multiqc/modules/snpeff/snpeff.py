#!/usr/bin/env python

""" MultiQC module to parse logs from SnpEff """

from __future__ import print_function

from collections import OrderedDict
import logging
from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ SnpEff """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SnpEff', anchor='snpeff',
            href="http://snpeff.sourceforge.net/",
            info=" is a genetic variant annotation and effect prediction toolbox. " \
                 "It annotates and predicts the effects of variants on genes (such as amino acid changes). ")

        self.snpeff_data = dict()
        self.snpeff_section_totals = dict()
        self.snpeff_qualities = dict()

        for f in self.find_log_files('snpeff', filehandles=True):
            self.parse_snpeff_log(f)

        # Filter to strip out ignored sample names
        self.snpeff_data = self.ignore_samples(self.snpeff_data)

        if len(self.snpeff_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.snpeff_data)))

        # Write parsed report data to a file
        self.write_data_file(self.snpeff_data, 'multiqc_snpeff')

        # General stats table
        self.general_stats()

        # Report sections
        self.add_section (
            name = 'Variants by Genomic Region',
            anchor = 'snpeff-genomic-regions',
            plot = self.count_genomic_region_plot()
        )
        self.add_section (
            name = 'Variant Effects by Impact',
            anchor = 'snpeff-effects-impact',
            plot = self.effects_impact_plot()
        )
        self.add_section (
            name = 'Variant Effects by Class',
            anchor = 'snpeff-functional-class',
            plot = self.effects_function_plot()
        )
        if len(self.snpeff_qualities) > 0:
            self.add_section (
                name = 'Variant Qualities',
                anchor = 'snpeff-qualities',
                plot = self.qualities_plot()
            )


    def parse_snpeff_log(self, f):
        """ Go through log file looking for snpeff output """

        keys = {
            '# Summary table': [
                'Genome', 'Number_of_variants_before_filter', 'Number_of_known_variants',
                'Number_of_effects', 'Genome_total_length', 'Change_rate'
            ],
            '# Effects by impact': [ 'HIGH', 'LOW', 'MODERATE', 'MODIFIER' ],
            '# Effects by functional class': [ 'MISSENSE', 'NONSENSE', 'SILENT', 'Missense_Silent_ratio' ],
            '# Hom/Het table': ['Het', 'Hom', 'Missing'],
            '# Ts/Tv summary': [ 'Transitions', 'Transversions', 'Ts_Tv_ratio' ],
            '# Count by effects': 'all',
            '# Count by genomic region': 'all'
        }
        parsed_data = {}
        section = None
        for l in f['f']:
            l = l.strip()
            if l[:1] == '#':
                section = l
                self.snpeff_section_totals[section] = dict()
                continue
            s = l.split(',')

            # Quality values / counts
            if section == '# Quality':
                quals = OrderedDict()
                if l.startswith('Values'):
                    values = [int(c) for c in l.split(',')[1:] ]
                    counts = f['f'].readline()
                    counts = [int(c) for c in counts.split(',')[1:] ]
                    c = 0
                    total = sum(counts)
                    for i, v in enumerate(values):
                        if c < (total * 0.995):
                            quals[v] = counts[i]
                            c += counts[i]
                if len(quals) > 0:
                    self.snpeff_qualities[f['s_name']] = quals

            # Everything else
            elif section in keys:
                if keys[section] == 'all' or any([k in s[0].strip() for k in keys[section]]):
                    try:
                        parsed_data[ s[0].strip() ] = float(s[1].strip())
                    except ValueError:
                        parsed_data[ s[0].strip() ] = s[1].strip()
                    except IndexError:
                        pass
                    else:
                        # Parsing the number worked - add to totals
                        try:
                            self.snpeff_section_totals[section][s[0].strip()] += parsed_data[ s[0].strip() ]
                        except KeyError:
                            self.snpeff_section_totals[section][s[0].strip()] = parsed_data[ s[0].strip() ]
                    if len(s) > 2 and s[2][-1:] == '%':
                        parsed_data[ '{}_percent'.format(s[0].strip()) ] = float(s[2][:-1])

        if len(parsed_data) > 0:
            if f['s_name'] in self.snpeff_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.snpeff_data[f['s_name']] = parsed_data

    def general_stats(self):
        """ Add key SnpEff stats to the general stats table """

        headers = OrderedDict()
        headers['Change_rate'] = {
            'title': 'Change rate',
            'scale': 'RdYlBu-rev',
            'min': 0,
            'format': '{:,.0f}'
        }
        headers['Ts_Tv_ratio'] = {
            'title': 'Ts/Tv',
            'description': 'Transitions / Transversions ratio',
            'format': '{:,.3f}'
        }
        headers['Number_of_variants_before_filter'] = {
            'title': 'M Variants',
            'description': 'Number of variants before filter (millions)',
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'min': 0,
            'format': '{:,.2f}'
        }
        self.general_stats_addcols(self.snpeff_data, headers)


    def count_genomic_region_plot(self):
        """ Generate the SnpEff Counts by Genomic Region plot """

        # Sort the keys based on the total counts
        keys = self.snpeff_section_totals['# Count by genomic region']
        sorted_keys = sorted(keys, reverse=True, key=keys.get)

        # Make nicer label names
        pkeys = OrderedDict()
        for k in sorted_keys:
            pkeys[k] = {'name': k.replace('_', ' ').title().replace('Utr', 'UTR') }

        # Config for the plot
        pconfig = {
            'id': 'snpeff_variant_effects_region',
            'title': 'SnpEff: Counts by Genomic Region',
            'ylab': '# Reads',
            'logswitch': True
        }

        return bargraph.plot(self.snpeff_data, pkeys, pconfig)


    def effects_impact_plot(self):
        """ Generate the SnpEff Counts by Effects Impact plot """

        # Put keys in a more logical order
        keys = [ 'MODIFIER', 'LOW', 'MODERATE', 'HIGH' ]

        # Make nicer label names
        pkeys = OrderedDict()
        for k in keys:
            pkeys[k] = {'name': k.title() }

        # Config for the plot
        pconfig = {
            'id': 'snpeff_variant_effects_impact',
            'title': 'SnpEff: Counts by Effects Impact',
            'ylab': '# Reads',
            'logswitch': True
        }

        return bargraph.plot(self.snpeff_data, pkeys, pconfig)

    def effects_function_plot(self):
        """ Generate the SnpEff Counts by Functional Class plot """

        # Cats to plot in a sensible order
        keys = [ 'SILENT', 'MISSENSE', 'NONSENSE' ]

        # Make nicer label names
        pkeys = OrderedDict()
        for k in keys:
            pkeys[k] = {'name': k.title() }

        # Config for the plot
        pconfig = {
            'id': 'snpeff_variant_effects_class',
            'title': 'SnpEff: Counts by Functional Class',
            'ylab': '# Reads',
            'logswitch': True
        }

        return bargraph.plot(self.snpeff_data, pkeys, pconfig)



    def qualities_plot(self):
        """ Generate the qualities plot """

        pconfig = {
            'smooth_points': 200,
            'id': 'snpeff_qualities',
            'title': 'SnpEff: Qualities',
            'ylab': 'Count',
            'xlab': 'Values',
            'xDecimals': False,
            'ymin': 0
        }

        return linegraph.plot(self.snpeff_qualities, pconfig)

