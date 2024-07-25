import logging
import os
from pathlib import Path
from typing import Optional

from multiqc.utils.util_functions import rmtree_with_retries

logger = logging.getLogger(__name__)

_tmp_dir: Optional[Path] = None


def get_tmp_dir() -> Path:
    """
    Lazy get or create tmp dir.
    Delays import of `tempfile` to allow setting `os.environ["TMPDIR"]` in tests
    from the `tmp_path` fixture, after `multiqc.report` is imported globally.
    """
    import tempfile

    global _tmp_dir
    if _tmp_dir is None:
        _tmp_dir = Path(tempfile.mkdtemp())
        logger.debug(f"Using new temporary directory: {_tmp_dir}")

    return _tmp_dir


def data_tmp_dir() -> Path:
    """
    Temporary directory to collect data files from running modules before copying to the final
    destination in multiqc.core.write_results
    """
    path = get_tmp_dir() / "multiqc_data"
    os.makedirs(path, exist_ok=True)
    return path


def plots_tmp_dir() -> Path:
    """
    Temporary directory to collect plot exports from running modules before copying to the final
    destination in multiqc.core.write_results
    """
    path = get_tmp_dir() / "multiqc_plots"
    os.makedirs(path, exist_ok=True)
    return path


def new_tmp_dir():
    global _tmp_dir
    _tmp_dir = None
