"""
Initialises when multiqc module is loaded.

Makes the following available under the main multiqc namespace:
- run()
- config
- __version__
"""

import sys
import warnings

warnings.filterwarnings("ignore", category=SyntaxWarning)

OLDEST_SUPPORTED_PYTHON_VERSION = "3.9"

if sys.version_info < tuple(map(int, OLDEST_SUPPORTED_PYTHON_VERSION.split("."))):
    raise RuntimeError(
        f"You are running MultiQC with Python {sys.version_info}. "
        f"Please upgrade Python! MultiQC does not support Python < {OLDEST_SUPPORTED_PYTHON_VERSION}, "
        "things will break."
    )

# Load config and report before anything else:
from multiqc import (
    config,
    report,
)
from multiqc.base_module import BaseMultiqcModule
from multiqc.interactive import (
    ClConfig,
    add_custom_content_section,
    get_general_stats_data,
    get_module_data,
    get_plot,
    list_data_sources,
    list_modules,
    list_plots,
    list_samples,
    load_config,
    parse_logs,
    reset,
    write_report,
)
from multiqc.multiqc import run
from multiqc.plots.plot import PConfig, Plot

__version__ = config.version

__all__ = [
    "run",
    "config",
    "report",
    "__version__",
    # The rest of the functions define the interactive use API:
    "parse_logs",
    "list_data_sources",
    "list_modules",
    "list_samples",
    "list_plots",
    "get_plot",
    "Plot",
    "PConfig",
    "get_module_data",
    "get_general_stats_data",
    "reset",
    "write_report",
    "add_custom_content_section",
    "BaseMultiqcModule",
    "load_config",
    "ClConfig",
]
