"""
Code to initialise the MultiQC logging
"""

import logging
import os
import shutil
import sys
import tempfile

import coloredlogs
import rich.console

from multiqc.utils import config, util_functions

log_tmp_dir = None
log_tmp_fn = "/dev/null"

rich_console: rich.console.Console


def is_running_in_notebook():
    try:
        from IPython import get_ipython

        if "IPKernelApp" in get_ipython().config:
            return True
    except (ImportError, AttributeError):
        pass
    return False


def init_log(logger, quiet: bool, verbose: int, no_ansi: bool = False):
    """
    Initializes logging.
    Prints logs to console with level defined by loglevel
    Also prints verbose log to the multiqc data directory if available.
    (multiqc_data/multiqc.log)

    loglevel (str): Determines the level of the log output.
    """
    log_level = "DEBUG" if verbose > 0 else "INFO"
    if quiet:
        log_level = "WARNING"
        config.quiet = True

    global rich_console
    rich_console = rich.console.Console(
        stderr=False,
        highlight=False,
        force_terminal=util_functions.force_term_colors(),
        color_system=None if no_ansi else "auto",
    )
    if not is_running_in_notebook():
        rich_console.print(
            f"\n  [dark_orange]///[/] [bold][link=https://multiqc.info]MultiQC[/link][/] :mag: [dim]| v{config.version}\n"
        )

    # File for logging
    global log_tmp_dir, log_tmp_fn
    # Have to create a separate directory for the log file otherwise Windows will complain
    # about same file used by different processes:
    log_tmp_dir = tempfile.mkdtemp()
    log_tmp_fn = os.path.join(log_tmp_dir, "multiqc.log")

    # Logging templates
    debug_template = "[%(asctime)s] %(name)-50s [%(levelname)-7s]  %(message)s"
    info_template = "|%(module)18s | %(message)s"

    # Remove log handlers left from previous calls to multiqc.run
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Base level setup
    logger.setLevel(getattr(logging, "DEBUG"))

    # Automatically set no_ansi if not a tty terminal
    if not no_ansi:
        if not sys.stderr.isatty() and not util_functions.force_term_colors():
            no_ansi = True

    # Set up the console logging stream
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(getattr(logging, log_level))
    level_styles = coloredlogs.DEFAULT_LEVEL_STYLES
    level_styles["debug"] = {"faint": True}
    field_styles = coloredlogs.DEFAULT_FIELD_STYLES
    field_styles["module"] = {"color": "blue"}
    if log_level == "DEBUG":
        if no_ansi:
            console.setFormatter(logging.Formatter(debug_template))
        else:
            console.setFormatter(
                coloredlogs.ColoredFormatter(fmt=debug_template, level_styles=level_styles, field_styles=field_styles)
            )
    else:
        if no_ansi:
            console.setFormatter(logging.Formatter(info_template))
        else:
            console.setFormatter(
                coloredlogs.ColoredFormatter(fmt=info_template, level_styles=level_styles, field_styles=field_styles)
            )
    logger.addHandler(console)

    # Now set up the file logging stream if we have a data directory
    file_handler = logging.FileHandler(log_tmp_fn, encoding="utf-8")
    file_handler.setLevel(getattr(logging, "DEBUG"))  # always DEBUG for the file
    file_handler.setFormatter(logging.Formatter(debug_template))
    logger.addHandler(file_handler)


def move_tmp_log():
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
