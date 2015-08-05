#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

import logging
import os
import shutil
import zipfile


class MultiqcModule(object):

    def __init__(self, analysis_dir, output_dir):

        # Static variables
        self.name = "FastQC"
        self.analysis_dir = analysis_dir
        self.output_dir = output_dir

        # Find and load any FastQC reports
        fastqc_data = {}
        for root, dirnames, filenames in os.walk(analysis_dir):
            # Extracted FastQC directory
            if root[-7:] == '_fastqc' and 'fastqc_data.txt' in filenames:
                s_name = os.path.basename(root)
                s_name = s_name[:-7]
                d_path = os.path.join(root, 'fastqc_data.txt')
                with open (d_path, "r") as f:
                    fastqc_data[s_name] = f.read()

            # Zipped FastQC report
            for f in filenames:
                if f[-11:] == '_fastqc.zip':
                    s_name = f[:-11]
                    d_name = f[:-4]
                    fqc_zip = zipfile.ZipFile(os.path.join(root, f))
                    with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as f:
                        fastqc_data[s_name] = f.read()

        if len(fastqc_data) == 0:
            logging.debug("Could not find any FastQC reports in {}".format(analysis_dir))
        else:
            logging.info("Found {} FastQC reports".format(len(fastqc_data)))
            # Copy across the required module files (CSS / Javascript etc)
            self.init_modfiles()

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
