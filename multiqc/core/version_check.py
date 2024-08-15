import logging
import os
import platform
import re
import sys
from typing import Optional, Dict, Union

import requests
from packaging import version

from multiqc import config
from multiqc.utils.util_functions import strtobool, is_running_in_notebook

logger = logging.getLogger(__name__)


def check_version(interactive_function_name: Optional[str] = None):
    # Check that we're running the latest version of MultiQC
    if config.no_version_check is True:
        return

    try:
        # Fetch the version info from the API
        meta: Dict[str, Union[str, bool, None]] = {
            "version_multiqc": config.short_version,
            "version_python": platform.python_version(),
            "operating_system": platform.system(),
            "is_docker": os.path.exists("/.dockerenv"),
            "is_singularity": os.path.exists("/.singularity.d"),
            "is_conda": os.path.exists(os.path.join(sys.prefix, "conda-meta")),
            "is_ci": strtobool(os.getenv("CI", False)),
            "is_notebook": is_running_in_notebook(),
            "interactive_function_name": interactive_function_name,
        }
        wait_seconds = 2
        try:
            r = requests.get(config.version_check_url, params=meta, timeout=(wait_seconds / 2, wait_seconds / 2))
        except requests.exceptions.Timeout as e:
            logger.debug(f"Timed out after waiting for {wait_seconds}s for multiqc.info to check latest version: {e}")
        except requests.exceptions.RequestException as e:
            logger.debug(f"Could not connect to multiqc.info for version check: {e}")
        else:
            release_info = r.json()
            # Broadcast log messages if found
            for msg in release_info.get("broadcast_messages", []):
                if msg.get("message"):
                    level = msg.get("level")
                    if level not in ["debug", "info", "warning", "error", "critical"]:
                        level = "info"
                    getattr(logger, level)(msg["message"])
            # Available update log if newer
            remove_version = version.parse(re.sub(r"[^0-9.]", "", release_info["latest_release"]["version"]))
            this_version = version.parse(re.sub(r"[^0-9.]", "", config.short_version))
            if remove_version > this_version:
                logger.warning(f"MultiQC Version {release_info['latest_release']['version']} now available!")
            logger.debug(
                f"Latest MultiQC version is {release_info['latest_release']['version']}, "
                f"released {release_info['latest_release']['release_date']}"
            )
    except Exception as e:
        logger.debug(f"Could not connect to multiqc.info for version check: {e}")
