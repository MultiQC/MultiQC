#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.

See the Click documentation for more command line flag types:
http://click.pocoo.org/5/
"""

import click

# Sets config.kwargs['disable_plugin'] to True if specified (will be False otherwise)

mgikit_disable_plugin = click.option(
    "--mgikit-disable-plugin",
    "mgikit_disable_plugin",
    is_flag=True,
    type=bool,
    default=False,
    help="Disable the mgikit MultiQC plugin on this run",
)

undetermined_barcode_threshold = click.option(
    "--mgikit-undetermined-barcode",
    "undetermined_barcode_threshold",
    type=int,
    default=25,
    help="Set the number of undetermined barcodes to be presented in the report.",
)


decimal_positions = click.option(
    "--mgikit-decimal-positions",
    "decimal_positions",
    type=int,
    default=2,
    help="The number of decimal positions to be used with the values in teh tables of mgikit. 2 decimal positions is the default option.",
)
brief_report = click.option(
    "--mgikit-brief-report",
    "brief_report",
    type=bool,
    default=False,
    is_flag=True,
    help="generate a brief version of the report. This ignores the reports for cluster per sample per lane.",
)

keep_core_samples = click.option(
    "--mgikit-keep-core-samples",
    "keep_core_samples",
    type=bool,
    default=False,
    is_flag=True,
    help="Ignore undetermined and ambiguous cases in the report.",
)
