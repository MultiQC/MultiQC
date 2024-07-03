import logging
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Optional

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
        logger.debug(f"Using temporary directory: {_tmp_dir}")

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


def rmtree_with_retries(path, _logger=None, max_retries=10):
    """
    Robustly tries to delete paths.
    Retries several times (with increasing delays) if an OSError
    occurs.  If the final attempt fails, the Exception is propagated
    to the caller.
    """

    for i in range(max_retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if _logger:
                _logger.info(f"Unable to remove path: {path}")
                _logger.info(f"Retrying after {i**2} seconds")
            else:
                print(f"Unable to remove path: {path}", file=sys.stderr)
                print(f"Retrying after {i**2} seconds", file=sys.stderr)
            time.sleep(i**2)

    # Final attempt, pass any Exceptions up to caller.
    shutil.rmtree(path)


def clean_up(_logger=None):
    global _tmp_dir
    if _tmp_dir is not None and Path(_tmp_dir).exists():
        rmtree_with_retries(_tmp_dir, _logger)
    _tmp_dir = None
