#!/usr/bin/env python3
"""
Fetch Material Design Icons SVG files from Iconify API.

This script fetches required SVG files from the Iconify API using the mdi (Material Design Icons)
icon set and saves them to the MultiQC source directory so they are bundled with the PyPI installation.
"""

import requests
from pathlib import Path
import logging
from typing import List

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Project root directory
ROOT_DIR = Path(__file__).parent.parent
MULTIQC_DIR = ROOT_DIR / "multiqc"
UTILS_DIR = MULTIQC_DIR / "utils"
ICONS_DIR = UTILS_DIR / "material_icons"

# Iconify API base URL
ICONIFY_API_BASE = "https://api.iconify.design"

# List of icons that MultiQC uses in mdi:iconname format
REQUIRED_ICONS: List[str] = [
    # Core functionality
    "mdi:information",
    "mdi:alert",
    "mdi:alert-circle",
    "mdi:help-circle",
    "mdi:clock",
    "mdi:magnify",
    # Navigation
    "mdi:chevron-down",
    "mdi:chevron-up",
    "mdi:chevron-left",
    "mdi:chevron-right",
    "mdi:arrow-left",
    "mdi:arrow-right",
    "mdi:menu",
    "mdi:close",
    # Toolbox actions
    "mdi:pin",
    "mdi:format-text",
    "mdi:eye",
    "mdi:eye-off",
    "mdi:download",
    "mdi:content-save",
    "mdi:school",
    "mdi:content-copy",
    # Table actions
    "mdi:table",
    "mdi:sort",
    "mdi:filter",
    "mdi:chart-scatter-plot",
    "mdi:equalizer",
    "mdi:format-align-left",
    "mdi:view-column",
    # File operations
    "mdi:folder",
    "mdi:file",
    "mdi:file-document",
    # Status indicators
    "mdi:check-circle",
    "mdi:cancel",
    "mdi:warning",
    "mdi:radiobox-blank",
    "mdi:radiobox-marked",
    "mdi:thumb-up",
    "mdi:thumb-down",
    # Color mode toggle icons
    "mdi:circle-half-full",
    "mdi:weather-sunny",
    "mdi:weather-night",
    "mdi:check",
    # Additional common icons
    "mdi:home",
    "mdi:view-dashboard",
    "mdi:cog",
    "mdi:refresh",
    "mdi:printer",
    "mdi:share",
    "mdi:link",
    "mdi:star",
    "mdi:bookmark",
    "mdi:heart",
    "mdi:pencil",
    "mdi:delete",
    "mdi:plus",
    "mdi:minus",
    "mdi:magnify-plus",
    "mdi:magnify-minus",
    "mdi:hand-pointing-up",
    # Additional icons found in current usage
    "mdi:image",
    "mdi:checkbox-marked",
    "mdi:checkbox-blank-outline",
    "mdi:folder-open",
    "mdi:content-save-outline",
    "mdi:qrcode",
    "mdi:format-bold",
    "mdi:violin",
    # About section icons
    "mdi:youtube",
    "mdi:book-open-page-variant",
    "mdi:github",
    "mdi:alert-circle-outline",
]


def fetch_icon_from_iconify(icon_name: str) -> bool:
    """
    Fetch an icon from Iconify API and save it as SVG.

    Args:
        icon_name: The iconify icon name (e.g., 'mdi:home')

    Returns:
        True if successful, False otherwise
    """
    # Extract icon set and name from iconify format (e.g., 'mdi:home' -> 'mdi', 'home')
    icon_set, icon_part = icon_name.split(":")
    url = f"{ICONIFY_API_BASE}/{icon_set}/{icon_part}.svg"

    # Convert to filename format (mdi:home -> mdi-home.svg)
    file_name = icon_name.replace(":", "-")
    dest_path = ICONS_DIR / f"{file_name}.svg"

    # Check if the file already exists
    if dest_path.exists():
        return True

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        # Save the SVG content
        with open(dest_path, "w", encoding="utf-8") as f:
            f.write(response.text)

        logger.info(f"Downloaded: {icon_name} -> {file_name}.svg")
        return True

    except requests.RequestException as e:
        logger.error(f"Failed to download {icon_name}: {e}")
        return False
    except Exception as e:
        logger.error(f"Error saving {icon_name}: {e}")
        return False


def fetch_icons():
    """Fetch Material Design Icons from Iconify API."""
    missing_icons = False

    # Create destination directory
    ICONS_DIR.mkdir(parents=True, exist_ok=True)

    # Fetch SVG files
    for icon_name in REQUIRED_ICONS:
        success = fetch_icon_from_iconify(icon_name)
        if not success:
            missing_icons = True

    # Create a simple LICENSE file for attribution
    license_content = """Material Design Icons

These icons are from the Material Design Icons collection.
Source: https://materialdesignicons.com/
License: Apache License 2.0 or SIL Open Font License 1.1
Downloaded via Iconify API: https://api.iconify.design/
"""

    license_path = ICONS_DIR / "LICENSE"
    with open(license_path, "w", encoding="utf-8") as f:
        f.write(license_content)

    return not missing_icons


if __name__ == "__main__":
    success = fetch_icons()
    if not success:
        exit(1)
    logger.info("Material Design Icons download completed successfully!")
