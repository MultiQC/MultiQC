"""
__init__.py
~~~~~~~~~~~~~~~~~~~~
Initialises when multiqc module is loaded.

Makes the following available under the main multiqc namespace:
- run()
- config
- config.logger
- __version__
"""

import logging

from .multiqc import (
    run,
    load,
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
from .utils import config

config.logger = logging.getLogger(__name__)

__version__ = config.version

__all__ = [
    "run",
    "load",  # for interactive use
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
