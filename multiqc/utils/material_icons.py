"""
Material Design Icons utility module for MultiQC.

This module provides functionality to load and use Material Design Icons
across all MultiQC templates and Python code using the Iconify naming scheme.
"""

import logging
from typing import Optional, Dict
from pathlib import Path
import re

logger = logging.getLogger(__name__)

# Cache for loaded SVG content
_svg_cache: Dict[str, str] = {}


def get_material_icon_path() -> Path:
    """Get the path to the Material Design Icons directory."""
    return Path(__file__).parent / "material_icons"


def load_svg_content(icon_name: str) -> Optional[str]:
    """
    Load SVG content for a given Material Design Icon.

    Args:
        icon_name: Name of the icon in mdi:name format (e.g., 'mdi:information', 'mdi:alert')

    Returns:
        SVG content as string, or None if not found
    """
    if icon_name in _svg_cache:
        return _svg_cache[icon_name]

    # Convert mdi:name to mdi-name for file path
    file_name = icon_name.replace(":", "-")
    icon_path = get_material_icon_path() / f"{file_name}.svg"

    if not icon_path.exists():
        logger.warning(f"Material Design Icon not found: {icon_name} (file: {file_name}.svg)")
        return None

    try:
        with open(icon_path, "r", encoding="utf-8") as f:
            svg_content = f.read()
            _svg_cache[icon_name] = svg_content
            return svg_content
    except Exception as e:
        logger.error(f"Failed to load Material Design Icon '{icon_name}': {e}")
        return None


def get_material_icon(
    icon_name: str, size: int = 24, color: Optional[str] = None, class_name: Optional[str] = None
) -> str:
    """
    Get a Material Design Icon as HTML SVG.

    Args:
        icon_name: Name of the icon in mdi:name format (e.g., 'mdi:information', 'mdi:alert')
        size: Size in pixels (default: 24)
        color: CSS color value (default: currentColor)
        class_name: Additional CSS class names

    Returns:
        HTML SVG string
    """
    svg_content = load_svg_content(icon_name)

    if svg_content is None:
        # Fallback to a basic icon if the requested icon is not found
        return f'<span class="material-icon-missing" title="Missing icon: {icon_name}">[{icon_name}]</span>'

    # Parse and modify the SVG
    # Replace width/height attributes
    svg_content = re.sub(r'width="[^"]*"', f'width="{size}"', svg_content)
    svg_content = re.sub(r'height="[^"]*"', f'height="{size}"', svg_content)

    # Add color and class attributes
    style_attrs = []
    if color:
        style_attrs.append(f'fill="{color}"')
    if class_name:
        # Add class attribute after the opening <svg tag
        svg_content = svg_content.replace("<svg", f'<svg class="{class_name}"', 1)

    # Add inline style if needed
    if style_attrs:
        style_str = " ".join(style_attrs)
        svg_content = svg_content.replace("<svg", f'<svg style="{style_str}"', 1)

    return svg_content


def get_material_icon_js(icon_name: str, size: int = 24, color: Optional[str] = None) -> str:
    """
    Get a Material Design Icon as JavaScript string for embedding in JS files.

    Args:
        icon_name: Name of the icon in mdi:name format
        size: Size in pixels (default: 24)
        color: CSS color value (default: currentColor)

    Returns:
        JavaScript string literal with SVG content
    """
    svg_content = get_material_icon(icon_name, size, color)
    # Escape for JavaScript
    js_content = svg_content.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
    return f'"{js_content}"'
