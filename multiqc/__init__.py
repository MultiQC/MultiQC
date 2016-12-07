import logging
import warnings

logger = logging.getLogger(__name__)

config = version = None

try:
    from pkg_resources import get_distribution
    from multiqc.utils import config
    __version__ = get_distribution("multiqc").version
    config.version = __version__
except:
    warnings.warn("Cannot determine MultiQC version. Configuration problem?")

