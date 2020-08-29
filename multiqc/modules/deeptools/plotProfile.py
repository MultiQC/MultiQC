#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotProfile """

import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class plotProfileMixin():
    def parse_plotProfile(self):
        """Find plotProfile output"""
        self.deeptools_plotProfile = dict()
        for f in self.find_log_files('deeptools/plotProfile', filehandles=False):
            parsed_data, bin_labels, converted_bin_labels = self.parsePlotProfileData(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotProfile:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotProfile[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotProfile')

        self.deeptools_plotProfile = self.ignore_samples(self.deeptools_plotProfile)

        if len(self.deeptools_plotProfile) > 0:

            # Try to do plot bands but don't crash if the labels aren't as we expect
            xPlotBands = []
            xPlotLines = []
            plotBandHelp = ''
            try:
                xPlotBands.append({'from': converted_bin_labels[bin_labels.index('TES')], 'to': converted_bin_labels[-1], 'color': '#f7cfcf'})
                xPlotBands.append({'from': converted_bin_labels[bin_labels.index('TSS')], 'to': converted_bin_labels[bin_labels.index('TES')], 'color': '#ffffe2'})
                xPlotBands.append({'from': converted_bin_labels[0], 'to': converted_bin_labels[bin_labels.index('TSS')], 'color': '#e5fce0'})
                xPlotLines.append({'width': 1, 'value': converted_bin_labels[bin_labels.index('TES')], 'dashStyle': 'Dash', 'color': '#000000'})
                xPlotLines.append({'width': 1, 'value': converted_bin_labels[bin_labels.index('TSS')], 'dashStyle': 'Dash', 'color': '#000000'})
                plotBandHelp = """
                    * Green: {} upstream of gene to {}
                    * Yellow: {} to {}
                    * Pink: {} to {} downstream of gene
                    """.format(
                        list(filter(None,bin_labels))[0],
                        list(filter(None,bin_labels))[1],
                        list(filter(None,bin_labels))[1],
                        list(filter(None,bin_labels))[2],
                        list(filter(None,bin_labels))[2],
                        list(filter(None,bin_labels))[3]
                    )
            except ValueError:
                pass


            config = {
                'id': 'read_distribution_profile',
                'title': 'deeptools: Read Distribution Profile after Annotation',
                'ylab': 'Occurrence',
                'xlab': None,
                'smooth_points': 100,
                'xPlotBands': xPlotBands,
                'xPlotLines': xPlotLines
            }

            self.add_section (
                name = 'Read Distribution Profile after Annotation',
                anchor = 'read_distribution_profile_plot',
                description = """
                    Accumulated view of the distribution of sequence reads related to the closest annotated gene.
                    All annotated genes have been normalized to the same size.

                    {}""".format(plotBandHelp),
                plot=linegraph.plot(self.deeptools_plotProfile, config)
            )

        return len(self.deeptools_plotProfile)

    def parsePlotProfileData(self, f):
        d = dict()
        bin_labels = []
        bins = []
        for line in f['f'].splitlines():
            cols = line.rstrip().split("\t")
            if cols[0] == "bin labels":
                for col in cols[2:len(cols)]:
                    if col not in list(filter(None,bin_labels)):
                        bin_labels.append(col)
                    else:
                        break
            elif cols[0] == "bins":
                for col in cols[2:len(cols)]:
                    if len(bins)!=len(bin_labels):
                        bins.append(self._int(col))
                    else:
                        break
            else:
                s_name = self.clean_s_name(cols[0], f['root'])
                d[s_name] = dict()

                factors = {'Kb': 1e3, 'Mb': 1e6, 'Gb': 1e9}
                convert_factor = 1
                for k, v in factors.items():
                    if k in bin_labels[0]:
                        convert_factor *= v
                        start = float(bin_labels[0].strip(k)) * convert_factor
                step = self._int(abs(start/bin_labels.index('TSS')))
                end = step*(len(bin_labels)-bin_labels.index('TSS')-1)
                converted_bin_labels = range((self._int(start)+step), (self._int(end)+step), step)

                for i in bins:
                    d[s_name].update({converted_bin_labels[i-1]:float(cols[i+1])})

        return d, bin_labels, converted_bin_labels
