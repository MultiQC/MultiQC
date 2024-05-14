"""
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

# Load config and report before anything else:
from multiqc import config  # noqa: E402
from multiqc import report  # noqa: E402
from multiqc.multiqc import run  # noqa: E402
from multiqc.base_module import BaseMultiqcModule  # noqa: E402
from multiqc.interactive import (  # noqa: E402
    parse_logs,
    list_data_sources,
    list_modules,
    list_samples,
    list_plots,
    get_plot,
    get_module_data,
    get_general_stats_data,
    reset,
    write_report,
    add_custom_content_section,
    load_config,
)
from multiqc.plots.plotly.plot import Plot, PConfig  # noqa: E402

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
]
