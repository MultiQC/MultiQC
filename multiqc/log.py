import os
import sys
import logging

LEVELS = {
    0 : 'WARNING',
    1 : 'INFO',
    2 : 'DEBUG',
}

def init_log(logger, filename=None, loglevel=None):
    """
    Initializes the log file in the proper format.

    Arguments:

        filename (str): Path to a file. Or None if logging is to
                         be disabled.
        loglevel (str): Determines the level of the log output.
    """
    template = '[%(asctime)s] %(levelname)-8s: %(name)-25s: %(message)s'
    formatter = logging.Formatter(template)

    if loglevel:
        logger.setLevel(getattr(logging, loglevel))

    # We will always print warnings and higher to stderr
    console = logging.StreamHandler()
    console.setLevel('WARNING')
    console.setFormatter(formatter)

    if filename:
        file_handler = logging.FileHandler(filename, encoding='utf-8')
        if loglevel:
            file_handler.setLevel(getattr(logging, loglevel))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # If no logfile is provided we print all log messages that the user has
    # defined to stderr
    else:
        if loglevel:
            console.setLevel(getattr(logging, loglevel))

    logger.addHandler(console)

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
