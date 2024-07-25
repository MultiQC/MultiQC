"""
Code to initialise the MultiQC logging
"""

import datetime
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import List, Optional

import coloredlogs  # type: ignore
import rich
from rich.logging import RichHandler
from rich.theme import Theme
import rich.jupyter
import rich_click
import rich.progress
from tqdm import tqdm

from multiqc import config
from multiqc.core import tmp_dir
from multiqc.utils import util_functions
from multiqc.utils.util_functions import is_running_in_notebook

log_tmp_fn: Optional[Path] = None
log_file_handler: Optional[logging.FileHandler] = None

rich_console: Optional[rich.console.Console] = None

logger = logging.getLogger()  # root logger


DEBUG_TEMPLATE = "[%(asctime)s] %(name)-50s [%(levelname)-7s]  %(message)s"


def init_log(log_to_file: bool = False):
    """
    Initializes logging.
    Prints logs to console with level defined by loglevel
    Also prints verbose log to the multiqc data directory if available.
    (multiqc_data/multiqc.log)

    loglevel (str): Determines the level of the log output.
    """
    # Remove log handlers left from previous calls to multiqc.run. Makes the function idempotent
    logger.handlers.clear()

    # Console log level
    log_level = "DEBUG" if config.verbose else "INFO"
    if config.quiet:
        log_level = "WARNING"
    logger.setLevel(log_level)

    # Remove DEBUG level for the PIL.PngImagePlugin logger
    logging.getLogger("PIL").setLevel(logging.INFO)

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
        stderr=True,
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

    _setup_coloredlogs(log_level, logger)

    if not config.quiet:
        if util_functions.is_running_in_notebook():
            _print_intro_with_coloredlogs()
        else:
            _print_intro_with_rich()

    if log_to_file:
        _add_file_handler()


def _add_file_handler() -> logging.Handler:
    """
    Set up the file logging stream if we have a data directory
    """
    global log_tmp_fn
    if log_tmp_fn is None:
        log_tmp_fn = tmp_dir.get_tmp_dir() / "multiqc.log"

    global log_file_handler
    log_file_handler = logging.FileHandler(log_tmp_fn, encoding="utf-8")
    log_file_handler.setLevel(logging.DEBUG)  # always DEBUG for the file
    log_file_handler.setFormatter(logging.Formatter(DEBUG_TEMPLATE))
    logger.addHandler(log_file_handler)
    logger.debug(f"Logging to file: {log_tmp_fn}")
    return log_file_handler


def _setup_coloredlogs(log_level, logger):
    # Use coloredlogs as Rich is breaking output formatting
    info_template = "%(module)18s | %(message)s"

    # Set up the console logging stream
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(log_level)
    level_styles = coloredlogs.DEFAULT_LEVEL_STYLES
    level_styles["debug"] = {"faint": True}
    field_styles = coloredlogs.DEFAULT_FIELD_STYLES
    field_styles["module"] = {"color": "blue"}
    field_styles["name"] = {"color": "blue"}

    if log_level == "DEBUG":
        if config.no_ansi is True:
            console.setFormatter(logging.Formatter(DEBUG_TEMPLATE))
        else:
            console.setFormatter(
                coloredlogs.ColoredFormatter(fmt=DEBUG_TEMPLATE, level_styles=level_styles, field_styles=field_styles)
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


def _setup_rich_handler(log_level, console, logger):
    # Set up the console logging stream
    console_handler = RichHandler(
        level=log_level,
        console=console,
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
    emoji = choose_emoji()
    emoji = f" {emoji}" if emoji else ""
    intro = f"{DARK_ORANGE}///{RESET} {BOLD}https://multiqc.info{RESET}{emoji} {DIM}v{config.version}{RESET}"
    if not util_functions.is_running_in_notebook():
        intro = f"\n{intro}\n"
    print(intro)


def _print_intro_with_rich():
    if rich_console is not None:
        rich_console.print(f"\n{rich_click.rich_click.HEADER_TEXT}\n")


def remove_file_handler():
    """
    Move the temporary log file to the MultiQC data directory if it exists, and remove the file handler.
    """

    # https://stackoverflow.com/questions/15435652/python-does-not-release-filehandles-to-logfile
    if log_file_handler is not None:
        log_file_handler.close()
        logger.removeHandler(log_file_handler)

    global log_tmp_fn
    if log_tmp_fn is not None:
        if log_tmp_fn.exists():
            if config.data_dir is not None and Path(config.data_dir).is_dir():
                try:
                    shutil.copy(log_tmp_fn, Path(config.data_dir) / "multiqc.log")
                except IOError:
                    pass
        try:
            os.remove(log_tmp_fn)
        except OSError:
            pass
        log_tmp_fn = None


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


emoji_rich_ids = {
    "🍾": ":bottle_with_popping_cork:",
    "🌹": ":rose:",
    "🍀": ":four_leaf_clover:",
    "🌏": ":globe_showing_asia-australia:",
    "🎃": ":jack-o-lantern:",
    "🎅": ":santa:",
    "🎄": ":christmas_tree:",
    "🔍": ":mag:",
}

emoji_dates = {
    "🍾": (1, 1, 1, 5),  # New Year's Day
    "🌹": (2, 14, 0, 0),  # Valentine's Day
    "🍀": (3, 17, 0, 0),  # St Patrick's Day
    "🌏": (4, 22, 0, 0),  # Earth Day
    "🎃": (10, 31, 5, 0),  # Halloween
    "🎅": (12, 25, 0, 0),  # Christmas Day
    "🎄": (12, 25, 7, 7),  # Christmas
}


def choose_emoji(use_rich=False) -> str:
    """Choose an emoji to use in the report header."""
    # NB: We haven't parsed the config yet, so can't disable via config
    if _no_unicode():
        return ""

    today = datetime.date.today()

    selected_emoji = "🔍"
    for emoji, (month, day, days_before, days_after) in emoji_dates.items():
        special_date = datetime.date(today.year, month, day)
        date_range_start = special_date - datetime.timedelta(days=days_before)
        date_range_end = special_date + datetime.timedelta(days=days_after)
        if date_range_start <= today <= date_range_end:
            selected_emoji = emoji
    if use_rich:
        return emoji_rich_ids[selected_emoji]
    return selected_emoji


def _no_unicode() -> bool:
    # When LANG or PYTHONIOENCODING or is not set, Rich won't be able to print fancy unicode
    # characters for the progress bar, and the runtime would crash with UnicodeEncodeError:
    # https://github.com/MultiQC/MultiQC/actions/runs/8814275065/job/24193771822
    # See https://github.com/Textualize/rich/issues/212
    return (
        "utf".casefold() not in os.environ.get("LANG", "").casefold()
        and "utf".casefold() not in os.environ.get("PYTHONIOENCODING", "").casefold()
    )


def iterate_using_progress_bar(items: List, desc: str, update_fn, item_to_str_fn=str, disable_progress=False):
    # GitHub actions doesn't understand ansi control codes to move the cursor,
    # so it prints each update ona a new line. Better disable it for CI.
    disable_progress = disable_progress or config.no_ansi or config.quiet or os.getenv("CI") is not None

    # Rich widgets do not look good in Jupyter, of it there is no unicode support.
    # Additionally, falling back to tqdm if rich_console was not initialized. That
    # happens when init_log.init_log() wasn't run, i.e in unit tests.
    if is_running_in_notebook() or _no_unicode() or rich_console is None:
        # ANSI escape code for dim text
        if not config.no_ansi:
            DIM = "\033[2m"
            BLUE = "\033[34m"
            RESET = "\033[0m"
        else:
            DIM = ""
            BLUE = ""
            RESET = ""

        bar_format = f"{BLUE}{desc:>17} {RESET}| " + "{bar:40} {percentage:3.0f}% {n_fmt}/{total_fmt} {desc}"

        # Set up the tqdm progress bar
        with tqdm(
            total=len(items),
            desc=desc,
            unit="file",
            file=sys.stdout,
            disable=disable_progress,
            bar_format=bar_format,
        ) as pbar:
            for i, item in enumerate(items):
                pbar.update(1)
                # Update the progress bar description with the file being searched
                pbar.set_description_str(f"{DIM}{item_to_str_fn(item)[-50:]}{RESET}")
                update_fn(i, item)

            # Clear the description after the loop is complete
            pbar.set_description("")
            pbar.refresh()
    else:
        N_SPACES_BEFORE_PIPE = 15  # to align bar desc with other log entries
        need_to_add_spaces = max(0, N_SPACES_BEFORE_PIPE - len(desc))
        progress_obj = rich.progress.Progress(
            "[blue][/]" + " " * need_to_add_spaces,
            rich.progress.SpinnerColumn(),
            "[blue]{task.description}[/] |",
            rich.progress.BarColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            "[green]{task.completed}/{task.total}",
            "[dim]{task.fields[s_fn]}[/]",
            console=rich_console,
            disable=disable_progress,
        )
        with progress_obj as progress:
            mqc_task = progress.add_task(desc, total=len(items), s_fn="")
            for i, item in enumerate(items):
                progress.update(mqc_task, advance=1, s_fn=item_to_str_fn(item)[-50:])
                update_fn(i, item)
            progress.update(mqc_task, s_fn="")


def rich_console_print(*args, **kwargs):
    if rich_console is not None:
        rich_console.print(*args, **kwargs)
    else:
        print(*args, **kwargs)
