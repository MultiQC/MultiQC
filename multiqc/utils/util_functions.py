"""MultiQC Utility functions, used in a variety of places."""

import logging
from collections import defaultdict, OrderedDict

import os
import shutil
import sys
import time
import datetime


logger = logging.getLogger(__name__)


def robust_rmtree(path, logger=None, max_retries=10):
    """Robustly tries to delete paths.
    Retries several times (with increasing delays) if an OSError
    occurs.  If the final attempt fails, the Exception is propagated
    to the caller.
    """

    for i in range(max_retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if logger:
                logger.info(f"Unable to remove path: {path}")
                logger.info(f"Retrying after {i**2} seconds")
            else:
                print(f"Unable to remove path: {path}", file=sys.stderr)
                print(f"Retrying after {i**2} seconds", file=sys.stderr)
            time.sleep(i**2)

    # Final attempt, pass any Exceptions up to caller.
    shutil.rmtree(path)


def strtobool(val) -> bool:
    """
    Replaces deprecated https://docs.python.org/3.9/distutils/apiref.html#distutils.util.strtobool
    The deprecation recommendation is to re-implement the function https://peps.python.org/pep-0632/

    ------------------------------------------------------------

    Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val_str = str(val).lower()
    if val_str in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val_str in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")


def choose_emoji():
    """Choose an emoji to use in the report header."""
    # NB: We haven't parsed the config yet, so can't disable via config
    today = datetime.date.today()
    emojis = {
        "bottle_with_popping_cork": (1, 1, 1, 5),  # New Year's Day
        "rose": (2, 14, 0, 0),  # Valentine's Day
        "four_leaf_clover": (3, 17, 0, 0),  # St Patrick's Day
        "globe_showing_asia-australia": (4, 22, 0, 0),  # Earth Day
        "jack-o-lantern": (10, 31, 5, 0),  # Halloween
        "santa": (12, 25, 0, 0),  # Christmas Day
        "christmas_tree": (12, 25, 7, 7),  # Christmas
    }
    for emoji, (month, day, days_before, days_after) in emojis.items():
        special_date = datetime.date(today.year, month, day)
        date_range_start = special_date - datetime.timedelta(days=days_before)
        date_range_end = special_date + datetime.timedelta(days=days_after)
        if date_range_start <= today <= date_range_end:
            return emoji
    return "mag"


def replace_defaultdicts(data):
    """
    Recursively replace dict-likes as dicts for nice yaml representation.
    """

    def _replace(obj):
        if isinstance(obj, (defaultdict, OrderedDict, dict)):
            return {k: _replace(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_replace(v) for v in obj]
        elif isinstance(obj, set):
            return {_replace(v) for v in obj}
        return obj

    return _replace(data)


def is_running_in_notebook() -> bool:
    try:
        from IPython import get_ipython

        if "IPKernelApp" in get_ipython().config:
            return True
    except (ImportError, AttributeError):
        pass
    return False


def no_unicode() -> bool:
    # When LANG or PYTHONIOENCODING or is not set, Rich won't be able to print fancy unicode
    # characters for the progress bar, and the runtime would crash with UnicodeEncodeError:
    # https://github.com/MultiQC/MultiQC/actions/runs/8814275065/job/24193771822
    # See https://github.com/Textualize/rich/issues/212
    return (
        "utf".casefold() not in os.environ.get("LANG", "").casefold()
        and "utf".casefold() not in os.environ.get("PYTHONIOENCODING", "").casefold()
    )
