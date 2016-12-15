#!/usr/bin/env python

"""
Code to initilise the MultiQC logging
"""

import logging
import os
import shutil
import tempfile

from multiqc.utils import config, util_functions

LEVELS = {0: 'INFO', 1: 'DEBUG'}
log_tmp_dir = tempfile.mkdtemp()
log_tmp_fn = os.path.join(log_tmp_dir, 'multiqc.log')

def init_log(logger, loglevel=0):
    """
    Initializes logging.
    Prints logs to console with level defined by loglevel
    Also prints verbose log to the multiqc data directory if available.
    (multiqc_data/.multiqc_log)

    Args:
        loglevel (str): Determines the level of the log output.
    """
    # Logging templates
    debug_template = '[%(asctime)s] %(name)-50s [%(levelname)-7s]  %(message)s'
    info_template = '[%(levelname)-7s] %(module)15s : %(message)s'

    # Base level setup
    logger.setLevel(getattr(logging, 'DEBUG'))

    # Set up the console logging stream
    console = logging.StreamHandler()
    console.setLevel(getattr(logging, loglevel))
    if loglevel == 'DEBUG':
        console.setFormatter(logging.Formatter(debug_template))
    else:
        console.setFormatter(logging.Formatter(info_template))
    logger.addHandler(console)

    # Now set up the file logging stream if we have a data directory
    file_handler = logging.FileHandler(log_tmp_fn, encoding='utf-8')
    file_handler.setLevel(getattr(logging, 'DEBUG')) # always DEBUG for the file
    file_handler.setFormatter(logging.Formatter(debug_template))
    logger.addHandler(file_handler)

def copy_tmp_log(logger):
    """ Copy the temporary log file to the MultiQC data directory
    if it exists. """

    try:
        shutil.copyfile(log_tmp_fn, os.path.join(config.data_dir, '.multiqc.log'))
        # https://stackoverflow.com/questions/15435652/python-does-not-release-filehandles-to-logfile
        logger.shutdown()
        util_functions.robust_rmtree(log_tmp_dir)
    except (AttributeError, TypeError, IOError):
        pass


def get_log_stream(logger):
    """
    Returns a stream to the root log file.
    If there is no logfile return the stderr log stream

    Returns:
        A stream to the root log file or stderr stream.
    """

    file_stream = None
    log_stream = None
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            file_stream = handler.stream
        else:
            log_stream = handler.stream

    if file_stream:
        return file_stream

    return log_stream
