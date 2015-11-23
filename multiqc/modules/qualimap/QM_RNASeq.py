#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Submodule to handle code for Qualimap BamQC """

from __future__ import print_function
from collections import OrderedDict
import io
import logging
import os

from collections import defaultdict

from multiqc import config, BaseMultiqcModule

def parse_reports(self):
    """ Find Qualimap RNASeq reports and parse their data """
    
    return None


def report_sections(self):
    """ Add results from Qualimap RNASeq parsing to the report """
    # Append to self.sections list
    
    return None



def stats_table(self):
    """ Take the parsed stats from the QualiMap RNASeq report and add them to the
    basic stats table at the top of the report """
    
    return None

