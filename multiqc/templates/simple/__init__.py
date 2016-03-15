"""
========
 simple
========

This theme attempts to generate a report which is a simple as possible:
that is, has as little JavaScript as possible.

"""
import os

from multiqc.utils import config

template_parent = 'default'

config.plots_force_flat = True

template_dir = os.path.dirname(__file__)
base_fn = 'base.html'
