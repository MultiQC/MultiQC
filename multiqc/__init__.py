"""
__init__.py
~~~~~~~~~~~~~~~~~~~~
Initialises when multiqc module is loaded.

Makes the following available under the main multiqc namespace:
- run()
- config
- __version__
"""

import sys

OLDEST_SUPPORTED_PYTHON_VERSION = "3.8"

if sys.version_info < tuple(map(int, OLDEST_SUPPORTED_PYTHON_VERSION.split("."))):
    raise RuntimeError(
        "You are running MultiQC with Python {}. "
        "Please upgrade Python! MultiQC does not support Python < {}, "
        "things will break.".format(sys.version_info, OLDEST_SUPPORTED_PYTHON_VERSION)
    )

from .multiqc import (  # noqa: E402
    run,
    parse_logs,
    list_modules,
    list_samples,
    list_plots,
    show_plot,
    get_module_data,
    get_general_stats_data,
    reset,
    write_report,
    add_custom_content_section,
)
from .utils import config  # noqa: E402

__version__ = config.version

__all__ = [
    "run",
    "parse_logs",  # for interactive use
    "list_modules",  # for interactive use
    "list_samples",  # for interactive use
    "list_plots",  # for interactive use
    "show_plot",  # for interactive use
    "get_module_data",  # for interactive use
    "get_general_stats_data",  # for interactive use
    "reset",  # for interactive use
    "write_report",  # for interactive use
    "add_custom_content_section",  # for interactive use
    "config",
    "__version__",
]
