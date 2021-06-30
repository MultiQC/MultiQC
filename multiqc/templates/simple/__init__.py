"""
========
 simple
========

This theme attempts to generate a report which is a simple
as possible - no JavaScript where possible. The resulting
report is hopefully suitable for e-mailing, converting to PDF
and printing.

"""
import os

from multiqc.utils import config

template_parent = "default"

config.plots_force_flat = True
config.simple_output = True

template_dir = os.path.dirname(__file__)
base_fn = "base.html"
