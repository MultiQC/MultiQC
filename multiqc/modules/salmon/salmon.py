#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function

import json
import logging
import os
from collections import OrderedDict

import numpy as np

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.salmon.gcmodel import GCModel
from multiqc.modules.salmon.seqmodel import SeqModel
from multiqc.plots import linegraph
from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
        href='http://combine-lab.github.io/salmon/',
        info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename(os.path.dirname(f['root']))
            s_name = self.clean_s_name(s_name, f['root'])
            self.salmon_meta[s_name] = json.loads(f['f'])
        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        self.salmon_gc = []
        self.salmon_seq3 = []
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename(os.path.dirname(f['root']))
                s_name = self.clean_s_name(s_name, f['root'])
                self.parse_gc_bias(f['root'])
                self.parse_seq_bias(f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed
        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning
        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
            'title': '% Aligned',
            'description': '% Mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['num_mapped'] = {
            'title': 'M Aligned',
            'description': 'Mapped reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Fragment Length Distribution',
            'ylab': 'Fraction',
            'xlab': 'Fragment Length (bp)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_fld, pconfig))
        self.plot_gc_bias()
        self.plot_seq_bias()

    def plot_gc_bias(self):
        pconfig = lambda x: {
            'smooth_points': 25,
            'id': 'salmon_gc_plot {}'.format(x),
            'title': 'GC Bias {}'.format(x),
            'ylab': 'Obs/Exp ratio',
            'xlab': 'bins',
            'ymin': 0,
            'xmin': 0,
        }
        qrt1, qrt2, qrt3 = {}, {}, {}
        exp_avg = np.zeros(shape=(3, 25))
        obs_avg = np.zeros(shape=(3, 25))
        low, medium, high = [], [], []
        sample_names = []
        complete_avgs = []
        for sample_name, sample_gc in self.salmon_gc:
            sample_names.append(sample_name)
            exp = np.multiply(np.array(sample_gc.exp_weights_)[:, np.newaxis], sample_gc.exp_)
            obs = np.multiply(np.array(sample_gc.obs_weights_)[:, np.newaxis], sample_gc.obs_)
            exp_avg += exp
            obs_avg += obs
            ratio = np.divide(obs, exp)
            low.append(ratio[0])
            medium.append(ratio[1])
            high.append(ratio[2])
            complete_avgs.append(np.average([ratio[0], ratio[1], ratio[2]], axis=1))
            qrt1[sample_name] = self.scale(ratio[0], 100)
            qrt2[sample_name] = self.scale(ratio[1], 100)
            qrt3[sample_name] = self.scale(ratio[2], 100)
        low_bias_coeff = np.corrcoef(low)
        medium_bias_coeff = np.corrcoef(medium)
        high_bias_coeff = np.corrcoef(high)
        complete_avgs_coeff = np.corrcoef(complete_avgs)
        ratio_avg = np.divide(obs_avg, exp_avg)
        low_bias = self.scale(ratio_avg[0], 100)
        med_bias = self.scale(ratio_avg[1], 100)
        high_bias = self.scale(ratio_avg[2], 100)
        avg_plot = {'low-bias': low_bias, 'medium-bias': med_bias, 'high-bias': high_bias}
        self.add_section(plot=linegraph.plot(qrt1, pconfig('Low')))
        self.add_section(plot=linegraph.plot(qrt2, pconfig('Medium')))
        self.add_section(plot=linegraph.plot(qrt3, pconfig('High')))
        self.add_section(plot=linegraph.plot(avg_plot, pconfig('Average')))

        self.add_section(plot=heatmap.plot(low_bias_coeff, sample_names))
        self.add_section(plot=heatmap.plot(medium_bias_coeff, sample_names))
        self.add_section(plot=heatmap.plot(high_bias_coeff, sample_names))
        self.add_section(plot=heatmap.plot(complete_avgs_coeff, sample_names))

    def parse_gc_bias(self, f_root):
        bias_dir = os.path.dirname(f_root)
        sample_name = os.path.basename(os.path.dirname(bias_dir))
        is_exp_gc_exists = os.path.exists(os.path.join(bias_dir, 'aux_info', 'exp_gc.gz'))
        is_obs_gc_exists = os.path.exists(os.path.join(bias_dir, 'aux_info', 'obs_gc.gz'))
        if is_exp_gc_exists and is_obs_gc_exists:
            gc = GCModel()
            gc.from_file(bias_dir)
            self.salmon_gc.append((sample_name, gc))

    def scale(self, ratios, fragment_len):
        scaling_factor = fragment_len / (len(ratios))
        scaled_result = {}
        for i, ratio in enumerate(ratios):
                scaled_result[i * scaling_factor] = ratio
        return scaled_result

    def seq_scale(self, ratios, fragment_len):
        scaling_factor = fragment_len / (len(ratios))
        scaled_result = {}
        midpoint = len(ratios)//2
        print(len(ratios)//2)
        j = -1*midpoint
        while(j < midpoint):
            scaled_result[j] = ratios[j+midpoint]
            j=j+1
        print(scaled_result)
        return scaled_result

    def parse_seq_bias(self,f_root):
        bias_dir = os.path.dirname(f_root)
        sample_name = os.path.basename(os.path.dirname(bias_dir))
        is_exp_3_exists = os.path.exists(os.path.join(bias_dir, 'aux_info', 'exp3_seq.gz'))
        is_obs_3_exists = os.path.exists(os.path.join(bias_dir, 'aux_info', 'obs3_seq.gz'))
        if is_exp_3_exists and is_obs_3_exists:
            seq = SeqModel()
            seq.from_file(bias_dir)
            self.salmon_seq3.append((sample_name,seq))
            print("Multiply")

    def plot_seq_bias(self):

            pconfig = lambda x: {
                'smooth_points': 25,
                'id': 'salmon_seq_plot {}'.format(x),
                'title': 'Sequence Bias {}'.format(x),
                'ylab': 'Obs/Exp ratio',
                'xlab': 'windows',
            }

            Hpconfig = lambda x: {
                'title': 'Sequence Bias {}'.format(x)
            }

            row3A = row5A = {}
            row3C = row5C = {}
            row3G = row5G = {}
            row3T = row5T = {}

            first3 = []
            second3 = []
            third3 = []
            fourth3 = []
            first5 = []
            second5 = []
            third5 = []
            fourth5 = []

            exp_avg3 = np.zeros(shape=(4, 9))
            obs_avg3 = np.zeros(shape=(4, 9))

            exp_avg5 = np.zeros(shape=(4, 9))
            obs_avg5 = np.zeros(shape=(4, 9))

            sample_names = []
            for sample_name, sample_seq in self.salmon_seq3:
                sample_names.append(sample_name)
                # ratio3 = np.divide(sample_seq.obs3_, sample_seq.exp3_)
                # ratio5 = np.divide(sample_seq.obs5_, sample_seq.exp5_)
                # result3 = {'A': self.scale(ratio3[0], 100), 'C': self.scale(ratio3[1], 100), 'G': self.scale(ratio3[2], 100), 'T': self.scale(ratio3[3], 100)}
                # result5 = {'A': self.scale(ratio5[0], 100), 'C': self.scale(ratio5[1], 100), 'G': self.scale(ratio5[2], 100), 'T': self.scale(ratio5[3], 100)}

                ratio3 = np.divide(sample_seq.obs3_, sample_seq.exp3_)
                ratio5 = np.divide(sample_seq.obs5_, sample_seq.exp5_)
                exp_avg3 += sample_seq.exp3_
                obs_avg3 += sample_seq.obs3_

                exp_avg5 += sample_seq.exp5_
                obs_avg5 += sample_seq.obs5_

                first3.append(ratio3[0])
                second3.append(ratio3[1])
                third3.append(ratio3[2])
                fourth3.append(ratio3[3])

                # complete_avgs.append(np.average([ratio[0], ratio[1], ratio[2]], axis=1))

                first5.append(ratio5[0])
                second5.append(ratio5[1])
                third5.append(ratio5[2])
                fourth5.append(ratio5[3])

                row3A[sample_name] = self.seq_scale(ratio3[0], 100)
                row5A[sample_name] = self.seq_scale(ratio5[0], 100)

                row3C[sample_name] = self.seq_scale(ratio3[1], 100)
                row5C[sample_name] = self.seq_scale(ratio5[1], 100)

                row3G[sample_name] = self.seq_scale(ratio3[2], 100)
                row5G[sample_name] = self.seq_scale(ratio5[2], 100)

                row3T[sample_name] = self.seq_scale(ratio3[3], 100)
                row5T[sample_name] = self.seq_scale(ratio5[3], 100)


            ratio_avg3 = np.divide(obs_avg3, exp_avg3)
            ratio_avg5 = np.divide(obs_avg5, exp_avg5)

            A3_coeff = np.corrcoef(first3)
            C3_coeff = np.corrcoef(second3)
            G3_coeff = np.corrcoef(third3)
            T3_coeff = np.corrcoef(fourth3)

            A5_coeff = np.corrcoef(first5)
            C5_coeff = np.corrcoef(second5)
            G5_coeff = np.corrcoef(third5)
            T5_coeff = np.corrcoef(fourth5)

            avg_plot3 = {'A': self.scale(ratio_avg3[0], 100), 'C': self.scale(ratio_avg3[1], 100), 'G': self.scale(ratio_avg3[2], 100), 'T': self.scale(ratio_avg3[3], 100) }
            avg_plot5 = {'A': self.scale(ratio_avg5[0], 100), 'C': self.scale(ratio_avg5[1], 100), 'G': self.scale(ratio_avg5[2], 100), 'T': self.scale(ratio_avg5[3], 100) }

            self.add_section(plot=linegraph.plot(row3A, pconfig('3A')))
            self.add_section(plot=linegraph.plot(row3C, pconfig('3C')))
            self.add_section(plot=linegraph.plot(row3G, pconfig('3G')))
            self.add_section(plot=linegraph.plot(row3T, pconfig('3T')))

            self.add_section(plot=linegraph.plot(avg_plot3, pconfig('Average 3')))

            #
            self.add_section(plot=linegraph.plot(row5A, pconfig('5A')))
            self.add_section(plot=linegraph.plot(row5C, pconfig('5C')))
            self.add_section(plot=linegraph.plot(row5G, pconfig('5G')))
            self.add_section(plot=linegraph.plot(row5T, pconfig('5T')))

            self.add_section(plot=linegraph.plot(avg_plot5, pconfig('Average 5')))


            self.add_section(plot=heatmap.plot(A3_coeff, sample_names,Hpconfig('HeatMapA')))
            self.add_section(plot=heatmap.plot(C3_coeff, sample_names))
            self.add_section(plot=heatmap.plot(G3_coeff, sample_names))
            self.add_section(plot=heatmap.plot(T3_coeff, sample_names))


            self.add_section(plot=heatmap.plot(A5_coeff, sample_names))
            self.add_section(plot=heatmap.plot(C5_coeff, sample_names))
            self.add_section(plot=heatmap.plot(G5_coeff, sample_names))
            self.add_section(plot=heatmap.plot(T5_coeff, sample_names))
