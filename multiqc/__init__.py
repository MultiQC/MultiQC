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

from .multiqc import run, load, list_modules, list_samples, list_plots, show_plot
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
    "config",
    "__version__",
]
