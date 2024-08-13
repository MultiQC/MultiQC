import warnings

# noinspection PyUnresolvedReferences
from multiqc.report import *  # noqa: F403

# Issue a deprecation warning
warnings.warn(
    "Importing 'report' from 'multiqc.utils' is deprecated and will be removed in a future release. "
    "Please use 'from multiqc import report' instead.",
    DeprecationWarning,
    stacklevel=2,
)
