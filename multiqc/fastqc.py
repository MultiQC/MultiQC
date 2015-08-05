#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""


class MultiqcModule(object):

    def __init__(self, output_dir):

        self.name = "FastQC"
