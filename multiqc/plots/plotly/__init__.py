"""
Check that Plotly version is supported, and show an error and exit if it isn't.
"""

import logging
import sys

from packaging import version

logger = logging.getLogger(__name__)


LATEST_SUPPORTED = "5.17"


def check_plotly_version():
    try:
        import plotly
    except ImportError:
        logger.error("ERROR: Could not import Plotly")
        sys.exit(1)
    if version.parse(plotly.__version__) < version.parse(LATEST_SUPPORTED):
        logger.error(
            f"ERROR: found Plotly version {plotly.__version__}, but MultiQC needs at least {LATEST_SUPPORTED}. "
            f'Please upgrade: pip install "plotly>={LATEST_SUPPORTED}"'
        )
        sys.exit(1)
