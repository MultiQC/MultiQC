"""MultiQC Utility functions, used in a variety of places."""

from typing import Dict

import array
import json
import logging
from collections import defaultdict, OrderedDict

import os
import shutil
import sys
import time
import datetime
import math

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


def choose_emoji(rich=False) -> str:
    """Choose an emoji to use in the report header."""
    # NB: We haven't parsed the config yet, so can't disable via config
    if no_unicode():
        return ""

    today = datetime.date.today()

    selected_emoji = "🔍"
    for emoji, (month, day, days_before, days_after) in emoji_dates.items():
        special_date = datetime.date(today.year, month, day)
        date_range_start = special_date - datetime.timedelta(days=days_before)
        date_range_end = special_date + datetime.timedelta(days=days_after)
        if date_range_start <= today <= date_range_end:
            selected_emoji = emoji
    if rich:
        return emoji_rich_ids[selected_emoji]
    return selected_emoji


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
        elif isinstance(obj, tuple):
            return tuple(_replace(v) for v in obj)
        return obj

    return _replace(data)


def dump_json(data, filehandle, **kwargs):
    """
    Recursively replace non-JSON-conforming NaNs and lambdas with None.
    Note that a custom JSONEncoder would not work for NaNs:
    https://stackoverflow.com/a/28640141
    """

    def replace_nan(obj):
        """
        Recursively replace NaNs and Infinities with None
        """
        # Do checking in order of likelihood of occurrence
        if isinstance(obj, float):
            if math.isnan(obj) or math.isinf(obj):
                return None
            return obj
        if isinstance(obj, (tuple, set)):
            # JSON only knows list so convert tuples and sets to list.
            obj = list(obj)
        if isinstance(obj, list):
            for i, item in enumerate(obj):
                if isinstance(item, float) and (math.isnan(item) or math.isinf(item)):
                    obj[i] = None
                elif isinstance(item, (dict, list, tuple, set)):
                    obj[i] = replace_nan(item)
            return obj
        if isinstance(obj, dict):
            for key, value in obj.items():
                if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
                    obj[key] = None
                elif isinstance(value, (dict, list, tuple, set)):
                    obj[key] = replace_nan(value)
            return obj
        return obj

    class JsonEncoderWithArraySupport(json.JSONEncoder):
        """
        Encode array.array instances to list. Use the default method
        for this as it gets called only when an array instance is encountered
        and is then immediately serialized into a string. This saves memory
        compared to unpacking all arrays to list at once.
        """

        def default(self, o):
            if isinstance(o, array.array):
                return replace_nan(o.tolist())
            if callable(o):
                return None
            return super().default(o)

    if filehandle:
        json.dump(replace_nan(data), filehandle, cls=JsonEncoderWithArraySupport, **kwargs)
    else:
        return json.dumps(replace_nan(data), cls=JsonEncoderWithArraySupport, **kwargs)


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


def compress_number_lists_for_json(obj):
    """
    Take an object that should be JSON and compress all the lists of integer
    and lists of float as array.array. This saves space and the arrays can
    easily be converted back again, using the dump_json function above.

    The technical explanation:
    A python list is an array of pointers to python objects:
    {
        list metadata including length
        a pointer to an array of pointers: [
            PyObject *
            PyObject *
            etc.
        ]
    }
    A python float is very simple and takes 24 bytes.
    {
        PyTypeObject *type
        Py_ssize_t refcount
        double the actual floating point.
    }
    A python integer is slightly more complicated, but similar to the float. It
    takes 28 bytes for 32-bit data, 32 bytes for 64-bit, 36 bytes for 96-bit etc.

    An array.array is more simple.
    {
        array metadata including length
        a pointer to an array of machine values: [
            double,
            double,
            double,
            etc.
        ]
        more metadata
    }
    Using 8-byte machine values rather than Python objects saves thus
    24 bytes per float.
    """
    if isinstance(obj, (list, tuple)):
        try:
            # Try integer list first, because it does not accept floats.
            return array.array("q", obj)
        except TypeError:
            pass
        try:
            return array.array("d", obj)
        except TypeError:
            return [compress_number_lists_for_json(v) for v in obj]
    if isinstance(obj, dict):
        return {k: compress_number_lists_for_json(v) for k, v in obj.items()}
    return obj


def update_dict(target: Dict, source: Dict, none_only=False):
    """Recursively updates nested dict d from nested dict u"""
    assert target is not None, source is not None
    for key, val in source.items():
        if isinstance(val, dict):
            target[key] = update_dict(target.get(key, {}), val)
        else:
            if not none_only or target.get(key) is None:
                if isinstance(val, list):
                    target[key] = val.copy()
                else:
                    target[key] = val
    return target
