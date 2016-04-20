#!/usr/bin/env python

""" MultiQC functions to plot a report beeswarm group """

import json
import logging
import os
import random

from multiqc.utils import report, config
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

