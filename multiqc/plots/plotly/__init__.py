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
        import plotly  # type: ignore
    except ImportError:
        logger.error("ERROR: Could not import Plotly")
        sys.exit(1)
    if version.parse(plotly.__version__) < version.parse(LATEST_SUPPORTED):
        logger.error(
            f"ERROR: found Plotly version {plotly.__version__}, but MultiQC needs at least {LATEST_SUPPORTED}. "
            f'Please upgrade: pip install "plotly>={LATEST_SUPPORTED}"'
        )
        sys.exit(1)


def determine_barplot_height(max_n_samples, max_bars_in_group=1, legend_height=None):
    """
    Used in bar and box plots
    """

    n_bars = max_n_samples
    n_bars *= max_bars_in_group

    def n_bars_to_size(n: int):
        if n >= 30:
            return 15
        if n >= 20:
            return 20
        if n >= 10:
            return 25
        if n >= 5:
            return 30
        return 35

    # Set height to be proportional to the number of samples
    height = n_bars * n_bars_to_size(n_bars)

    if legend_height is not None:
        # expand the plot to fit the legend:
        height = max(height, max(height, legend_height))
        # but not too much - if there are only 2 samples, we don't want the plot to be too high:
        height = min(660, max(height, legend_height))

    height += 140  # Add space for the title and footer

    # now, limit the max and min height (plotly will start to automatically skip
    # some of the ticks on the left when the plot is too high)
    MIN_HEIGHT = 300
    MAX_HEIGHT = 2560
    height = min(MAX_HEIGHT, height)
    height = max(MIN_HEIGHT, height)
    return height
