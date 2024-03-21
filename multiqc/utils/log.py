"""
Code to initialise the MultiQC logging
"""

import logging
import os
import shutil
import sys
import tempfile
import rich
from rich.logging import RichHandler
from rich.theme import Theme

from multiqc.utils import config, util_functions

log_tmp_dir = None
log_tmp_fn = "/dev/null"

rich_console: rich.console.Console


def init_log(logger, verbose: int, quiet: bool, no_ansi: bool = False):
    """
    Initializes logging.
    Prints logs to console with level defined by loglevel
    Also prints verbose log to the multiqc data directory if available.
    (multiqc_data/multiqc.log)

    loglevel (str): Determines the level of the log output.
    """
    # File for logging
    global log_tmp_dir, log_tmp_fn
    log_tmp_dir = tempfile.mkdtemp()
    log_tmp_fn = os.path.join(log_tmp_dir, "multiqc.log")

    # Remove log handlers left from previous calls to multiqc.run
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Base level setup: used both for console and file logging
    logger.setLevel(getattr(logging, "DEBUG"))

    # Console log level
    log_level = "DEBUG" if verbose > 0 else "INFO"
    if quiet:
        log_level = "WARNING"
        config.quiet = True

    # Automatically set no_ansi if not a tty terminal
    if not no_ansi:
        if not sys.stderr.isatty() and not util_functions.force_term_colors():
            no_ansi = True

    # Set up the rich console
    global rich_console
    rich_console = rich.console.Console(
        stderr=True,
        highlight=False,
        force_terminal=util_functions.force_term_colors(),
        color_system=None if no_ansi else "auto",
        theme=Theme(
            styles={
                "logging.level.info": "",
                "logging.level.debug": "dim",
                "logging.level.warning": "yellow",
                "logging.level.error": "red",
                "repr.number": "",
                "repr.none": "",
                "repr.str": "",
                "repr.bool_true": "",
                "repr.bool_false": "",
                "log.time": "green",
            }
        ),
    )

    # Set up the console logging stream
    console_handler = RichHandler(
        level=log_level,
        console=rich_console,
        show_level=log_level == "DEBUG",
        show_time=log_level == "DEBUG",
        log_time_format="[%Y-%m-%d %H:%M:%S]",
        markup=True,
        omit_repeated_times=False,
        show_path=False,
    )

    class DebugFormatter(logging.Formatter):
        def format(self, record):
            if record.levelno == logging.DEBUG:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  [logging.level.debug]%(message)s[/]")
            elif record.levelno == logging.WARNING:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  [logging.level.warning]%(message)s[/]")
            elif record.levelno == logging.ERROR:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  [logging.level.error]%(message)s[/]")
            else:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  %(message)s")
            return logging.Formatter.format(self, record)

    class InfoFormatter(logging.Formatter):
        def format(self, record):
            if record.levelno == logging.WARNING:
                self._style = logging.PercentStyle("[blue]|%(module)18s[/] | [logging.level.warning]%(message)s[/]")
            elif record.levelno == logging.ERROR:
                self._style = logging.PercentStyle("[blue]|%(module)18s[/] | [logging.level.error]%(message)s[/]")
            else:
                self._style = logging.PercentStyle("[blue]|%(module)18s[/] | %(message)s")
            return logging.Formatter.format(self, record)

    console_handler.setLevel(getattr(logging, log_level))
    if log_level == "DEBUG":
        console_handler.setFormatter(DebugFormatter())
    else:
        console_handler.setFormatter(InfoFormatter())
    logger.addHandler(console_handler)

    # Now set up the file logging stream if we have a data directory
    file_handler = logging.FileHandler(log_tmp_fn, encoding="utf-8")
    file_handler.setLevel(getattr(logging, "DEBUG"))  # always DEBUG for the file
    file_handler.setFormatter(logging.Formatter("[%(asctime)s] %(name)-50s [%(levelname)-7s]  %(message)s"))
    logger.addHandler(file_handler)


def move_tmp_log(logger):
    """Move the temporary log file to the MultiQC data directory
    if it exists."""

    try:
        # https://stackoverflow.com/questions/15435652/python-does-not-release-filehandles-to-logfile
        logging.shutdown()
        shutil.copy(log_tmp_fn, os.path.join(config.data_dir, "multiqc.log"))
        os.remove(log_tmp_fn)
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
