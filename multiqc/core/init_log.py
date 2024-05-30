"""
Code to initialise the MultiQC logging
"""

import logging
import os
import shutil
import sys
import tempfile

import coloredlogs
import rich
from rich.logging import RichHandler
from rich.theme import Theme
import rich.jupyter
import rich_click

from multiqc import config
from multiqc.utils import util_functions


log_tmp_dir = None
log_tmp_fn = "/dev/null"

rich_console: rich.console.Console


def init_log():
    """
    Initializes logging.
    Prints logs to console with level defined by loglevel
    Also prints verbose log to the multiqc data directory if available.
    (multiqc_data/multiqc.log)

    loglevel (str): Determines the level of the log output.
    """
    global log_tmp_dir, log_tmp_fn
    # Have to create a separate directory for the log file otherwise Windows will complain
    # about same file used by different processes:
    log_tmp_dir = tempfile.mkdtemp()
    log_tmp_fn = os.path.join(log_tmp_dir, "multiqc.log")

    logger = logging.getLogger()  # root logger

    # Remove log handlers left from previous calls to multiqc.run. Makes the function idempotent
    logger.handlers.clear()

    # Console log level
    log_level = "DEBUG" if config.verbose else "INFO"
    if config.quiet:
        log_level = "WARNING"
    logger.setLevel(log_level)

    # Automatically set no_ansi if not a tty terminal
    if config.no_ansi is False:
        if not sys.stderr.isatty() and not force_term_colors():
            config.no_ansi = True

    # Reset margin-bottom to remove the gian gap between lines.
    # See https://github.com/Textualize/rich/issues/3335 for more context
    rich.jupyter.JUPYTER_HTML_FORMAT = rich.jupyter.JUPYTER_HTML_FORMAT.replace('style="', 'style="margin-bottom:0px;')

    # Set up the rich console
    global rich_console
    rich_console = rich.console.Console(
        stderr=False,
        highlight=False,
        force_terminal=force_term_colors(),
        force_interactive=False if config.no_ansi else None,
        color_system=None if config.no_ansi is True else "auto",
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

    debug_template = "[%(asctime)s] %(name)-50s [%(levelname)-7s]  %(message)s"
    _setup_coloredlogs(log_level, logger, debug_template)

    if not config.quiet:
        if util_functions.is_running_in_notebook():
            _print_intro_with_coloredlogs()
        else:
            _print_intro_with_rich()

    # Now set up the file logging stream if we have a data directory
    file_handler = logging.FileHandler(log_tmp_fn, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)  # always DEBUG for the file
    file_handler.setFormatter(logging.Formatter(debug_template))
    logger.addHandler(file_handler)


def _setup_coloredlogs(log_level, logger, debug_template):
    # Use coloredlogs as Rich is breaking output formatting
    info_template = "%(module)18s | %(message)s"

    # Set up the console logging stream
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(log_level)
    level_styles = coloredlogs.DEFAULT_LEVEL_STYLES
    level_styles["debug"] = {"faint": True}
    field_styles = coloredlogs.DEFAULT_FIELD_STYLES
    field_styles["module"] = {"color": "blue"}
    field_styles["name"] = {"color": "blue"}

    if log_level == "DEBUG":
        if config.no_ansi is True:
            console.setFormatter(logging.Formatter(debug_template))
        else:
            console.setFormatter(
                coloredlogs.ColoredFormatter(fmt=debug_template, level_styles=level_styles, field_styles=field_styles)
            )
    else:
        if config.no_ansi is True:
            console.setFormatter(logging.Formatter(info_template))
        else:
            console.setFormatter(
                coloredlogs.ColoredFormatter(fmt=info_template, level_styles=level_styles, field_styles=field_styles)
            )
    logger.addHandler(console)

    # Google Colab notebooks duplicate log messages without this, see
    # https://stackoverflow.com/a/55877763/341474
    logger.propagate = False


def _setup_rich_handler(log_level, logger):
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
            elif record.levelno == logging.CRITICAL:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  [logging.level.critical]%(message)s[/]")
            else:
                self._style = logging.PercentStyle("[blue]%(name)-50s[/]  %(message)s")
            return logging.Formatter.format(self, record)

    class InfoFormatter(logging.Formatter):
        def format(self, record):
            if record.levelno == logging.WARNING:
                self._style = logging.PercentStyle("[blue]%(module)18s[/] | [logging.level.warning]%(message)s[/]")
            elif record.levelno == logging.ERROR:
                self._style = logging.PercentStyle("[blue]%(module)18s[/] | [logging.level.error]%(message)s[/]")
            else:
                self._style = logging.PercentStyle("[blue]%(module)18s[/] | %(message)s")
            return logging.Formatter.format(self, record)

    console_handler.setLevel(log_level)
    if log_level == "DEBUG":
        console_handler.setFormatter(DebugFormatter())
    else:
        console_handler.setFormatter(InfoFormatter())
    logger.addHandler(console_handler)


def _print_intro_with_coloredlogs():
    # Print intro
    if config.no_ansi is False:
        BOLD = "\033[1m"
        DIM = "\033[2m"
        DARK_ORANGE = "\033[38;5;208m"  # ANSI code for dark orange color
        RESET = "\033[0m"
    else:
        BOLD = ""
        DIM = ""
        DARK_ORANGE = ""
        RESET = ""
    emoji = util_functions.choose_emoji()
    emoji = f" {emoji}" if emoji else ""
    intro = f"{DARK_ORANGE}///{RESET} {BOLD}https://multiqc.info{RESET}{emoji} {DIM}v{config.version}{RESET}"
    if not util_functions.is_running_in_notebook():
        intro = f"\n{intro}\n"
    print(intro)


def _print_intro_with_rich():
    rich_console.print(f"\n{rich_click.rich_click.HEADER_TEXT}\n")


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


def force_term_colors():
    """
    Check if any environment variables are set to force Rich to use coloured output
    """
    if os.getenv("GITHUB_ACTIONS") or os.getenv("FORCE_COLOR") or os.getenv("PY_COLORS"):
        return True
    return None
