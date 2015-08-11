#!/usr/bin/env python

""" MultiQC module to parse output from FastQ Screen """

import json
import logging
import os
import re

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "FastQ Screen"
        self.anchor = "fastq_screen"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/" target="_blank">FastQ Screen</a> \
            allows you to screen a library of sequences in FastQ format against a set of sequence databases so you can see if the \
            composition of the library matches with what you expect.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any FastQ Screen reports
        fq_screen_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith("_screen.txt"):
                    s_name = fn[:-11]
                    with open (os.path.join(root,fn), "r") as f:
                        fq_screen_raw_data[s_name] = f.read()

        if len(fq_screen_raw_data) == 0:
            logging.debug("Could not find any FastQ Screen reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} FastQ Screen reports".format(len(fq_screen_raw_data)))

        self.sections = list()

        # Section 1 - Alignment Profiles
        length_trimmed = self.parse_fqscreen(fq_screen_raw_data)
        self.sections.append({
            'name': 'Trimming Length Profiles',
            'anchor': 'fq_screen-lengths',
            'content': self.trimgalore_length_trimmed_plot(length_trimmed)
        })


    def parse_fqscreen(self, fq_screen_raw_data):
        """ Parse the FastQ Screen output into a 3D dict """
        parsed_data = {}
        for s, data in fq_screen_raw_data.iteritems():
            parsed_data[s] = {}
            for l in data.splitlines():
                if l[:18] == "%Hit_no_libraries:":
                    parsed_data[s]['No hits'] = {'percentages':{}}
                    parsed_data[s]['No hits']['percentages']['one_hit_one_library'] = float(l[19:])
                else:
                    fqs = re.search(r"^(\w+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
                    if fqs:
                        org = fqs.group(1)
                        parsed_data[s][org] = {'percentages':{}, 'counts':{}}
                        parsed_data[s][org]['counts']['reads_processed'] = int(r_seqs.group(2))
                        parsed_data[s][org]['counts']['unmapped'] = int(r_seqs.group(3))
                        parsed_data[s][org]['percentages']['unmapped'] = float(r_seqs.group(4))
                        parsed_data[s][org]['counts']['one_hit_one_library'] = int(r_seqs.group(5))
                        parsed_data[s][org]['percentages']['one_hit_one_library'] = float(r_seqs.group(6))
                        parsed_data[s][org]['counts']['multiple_hits_one_library'] = int(r_seqs.group(7))
                        parsed_data[s][org]['percentages']['multiple_hits_one_library'] = float(r_seqs.group(8))
                        parsed_data[s][org]['counts']['one_hit_multiple_libraries'] = int(r_seqs.group(9))
                        parsed_data[s][org]['percentages']['one_hit_multiple_libraries'] = float(r_seqs.group(10))
                        parsed_data[s][org]['counts']['multiple_hits_multiple_libraries'] = int(r_seqs.group(11))
                        parsed_data[s][org]['percentages']['multiple_hits_multiple_libraries'] = float(r_seqs.group(12))
        return parsed_data

    def fqscreen_plot (self, parsed_data):

        categories: ['Spruce','Human','Human_chrX','Human_chrY','Mouse']

        {
            name: 'multiple_hits_multiple_libraries',
            data: [5, 3, 4, 7, 2],
            stack: 'sample_1'
        }

        data = list()
        names = ['one_hit_one_library', 'multiple_hits_one_library', 'one_hit_multiple_libraries', 'multiple_hits_multiple_libraries']
        for n in names:
            n_data = list()
            for s in sorted(parsed_data):
                for org in s:
                    cats.append(org)
                    for atype in parsed_data[s][org]['percentages']:
                        pairs.append([l, parsed_data[s][l]['obs_exp']])
                        data[atype].append({
                            'name': atype,
                            'data': pairs,
                            'stack': s
                        })
            data.append({
                'name': n,
                'data': n_data,
                'stack': s
            })

        html = '<div id="fq_screen_plot" style="height:500px;"></div>'

        return html
