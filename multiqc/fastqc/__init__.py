#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

import collections
import logging
import os
import re
import shutil
import zipfile

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, analysis_dir, output_dir):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "FastQC"
        self.analysis_dir = analysis_dir
        self.output_dir = output_dir

        # Find and load any FastQC reports
        fastqc_raw_data = {}
        for root, dirnames, filenames in os.walk(analysis_dir):
            # Extracted FastQC directory
            if root[-7:] == '_fastqc' and 'fastqc_data.txt' in filenames:
                s_name = os.path.basename(root)
                s_name = s_name[:-7]
                d_path = os.path.join(root, 'fastqc_data.txt')
                with open (d_path, "r") as f:
                    fastqc_raw_data[s_name] = f.read()

            # Zipped FastQC report
            for f in filenames:
                if f[-11:] == '_fastqc.zip':
                    s_name = f[:-11]
                    d_name = f[:-4]
                    fqc_zip = zipfile.ZipFile(os.path.join(root, f))
                    try:
                        with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                            fastqc_raw_data[s_name] = f.read()
                    except KeyError:
                        pass # Can't find fastqc_raw_data.txt in the zip file

        if len(fastqc_raw_data) == 0:
            logging.debug("Could not find any FastQC reports in {}".format(analysis_dir))
            return

        logging.info("Found {} FastQC reports".format(len(fastqc_raw_data)))

        self.sections = ['fubar!']

        # Section 1 - Basic Stats
        parsed_stats = self.fastqc_basic_stats(fastqc_raw_data)
        parsed_stats = collections.OrderedDict(sorted(parsed_stats.items()))
        stats_table_headers = {
            'percent_duplicates': '% Sequence Duplication',
            'sequence_length': 'Sequence Length (bp)',
            'percent_gc': '% GC',
            'total_sequences_m': 'Total Sequences (millions)'
        }
        self.sections.append({
            'name': 'Basic Stats',
            'content': self.dict_to_table(parsed_stats, colheaders=stats_table_headers, sort_rows=True)
        })


        # Copy across the required module files (CSS / Javascript etc)
        self.init_modfiles()

    def fastqc_basic_stats(self, fastqc_raw_data):
        """ Parse the contents of multiple fastqc_raw_data.txt files.
        Returns a 2D dict with basic stats, sample names as first keys,
        then statistic type as second key. """
        parsed_data = {}
        for s, data in fastqc_raw_data.iteritems():
            parsed_data[s] = {}

            dups = re.search("#Total Deduplicated Percentage\s+([\d\.]+)", data)
            if dups:
                parsed_data[s]['percent_duplicates'] = "{0:.2f}%".format(float(dups.group(1)))

            seqlen = re.search("Sequence length\s+(\d+)", data)
            if seqlen:
                parsed_data[s]['sequence_length'] = seqlen.group(1)

            gc = re.search("%GC\s+(\d+)", data)
            if gc:
                parsed_data[s]['percent_gc'] = "{}%".format(gc.group(1))

            numseq = re.search("Total Sequences\s+(\d+)", data)
            if numseq:
                parsed_data[s]['total_sequences_m'] = "{0:.1f}".format(float(numseq.group(1))/1000000)

        return parsed_data


    def init_modfiles(self):
        """ Copy the required assets into the output directory.
        Need to do this manually as we don't want to over-write
        the existing directory tree."""

        self.css = [ os.path.join('assets', 'css', 'multiqc_fastqc.css') ]
        self.js = [ os.path.join('assets', 'js', 'multiqc_fastqc.js') ]

        for f in self.css + self.js:
            d = os.path.join(self.output_dir, os.path.dirname(f))
            if not os.path.exists(d):
                os.makedirs(d)
            if not os.path.exists(os.path.join(self.output_dir, f)):
                shutil.copy(os.path.join(os.path.dirname(__file__), f), os.path.join(self.output_dir, f))
