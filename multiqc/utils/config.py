import warnings

# noinspection PyUnresolvedReferences
from multiqc.config import *  # noqa: F403

# Issue a deprecation warning
warnings.warn(
    "Importing 'config' from 'multiqc.utils' is deprecated and will be removed in a future release. "
    "Please use 'from multiqc import config' instead.",
    DeprecationWarning,
    stacklevel=2,
)
