#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" __main__.py
~~~~~~~~~~~~~~~~~~~~
Called when multiqc namespace is used, primarily by the cli:

$ multiqc .
$ python -m multiqc .
"""

from __future__ import print_function
import click
import pkg_resources
import multiqc


def modify_usage_error(main_command):
    ''' Function to modify the default click error handling.
    Used here to tell the user about how to find additional help.
    With thanks to this Stack Overflow answer: http://stackoverflow.com/a/43922088/713980
    :param main_command: top-level group or command object constructed by click wrapper
    :return: None
    '''
    def show(self, file=None):
        if file is None:
            file = click._compat.get_text_stderr()
        color = None
        if self.ctx is not None:
            color = self.ctx.color
            click.utils.echo(self.ctx.get_usage() + '\n', file=file, color=color)
        click.utils.echo('Error: %s\n\nThis is MultiQC v{}\n\nFor more help, run \'multiqc --help\' or visit http://multiqc.info\n'.format(multiqc.__version__) % self.format_message(), file=file, color=color)
    click.exceptions.UsageError.show = show


if __name__ == "__main__" or __name__ == 'multiqc.__main__':
    # Add any extra plugin command line options
    for entry_point in pkg_resources.iter_entry_points('multiqc.cli_options.v1'):
        opt_func = entry_point.load()
        multiqc = opt_func(multiqc)
    # Modify the default click error handling
    modify_usage_error(multiqc)
    # Call the main function
    multiqc.multiqc.run_cli(prog_name='multiqc')
