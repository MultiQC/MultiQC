#!/usr/bin/env python

"""
Called when multiqc namespace is used, primarily by the cli:

$ multiqc .
$ python -m multiqc .
"""

from importlib_metadata import entry_points

from . import multiqc


def run_multiqc():
    # Add any extra plugin command line options
    for entry_point in entry_points(group="multiqc.cli_options.v1"):
        opt_func = entry_point.load()
        multiqc.run_cli = opt_func(multiqc.run_cli)
    # Call the main function
    multiqc.run_cli(prog_name="multiqc")


# Script is run directly
# NB: Usually runs with pyproject.toml console_scripts: multiqc = "multiqc.__main__:run_multiqc"
if __name__ == "__main__":
    run_multiqc()
