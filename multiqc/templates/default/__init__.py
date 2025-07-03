"""
=========
 default
=========

The main MultiQC report template, lovingly known to its admirers
simply as "default"

Note, this is where most of the MultiQC report interactive functionality
is based and will be developed. Unless you want to do some really radical
changes, you probably don't want to replace this theme. Instead, you can
create a child theme that starts with 'default' and then overwrites
certain files.

For more information about creating child themes, see the docs:
docs/templates.md

"""

import os

from multiqc.core.strict_helpers import lint_error

template_dir = os.path.dirname(__file__)
base_fn = "base.html"


def get_material_icon(icon_name, size=24, color="currentColor"):
    """
    Load Material Design Icon SVG from node_modules.

    Args:
        icon_name: Name of the icon (e.g., 'delete', 'info', 'warning')
        size: Size of the icon in pixels (default: 24)
        color: Color of the icon (default: 'currentColor')

    Returns:
        SVG string or empty string if icon not found
    """

    def process_svg(svg_content):
        """Process SVG content with size, color, and accessibility attributes."""
        # Replace default attributes
        svg_content = svg_content.replace('width="24"', f'width="{size}"')
        svg_content = svg_content.replace('height="24"', f'height="{size}"')
        # Add color if not already present
        if "fill=" not in svg_content:
            svg_content = svg_content.replace("<svg", f'<svg fill="{color}"')
        # Add aria-hidden for accessibility
        svg_content = svg_content.replace("<svg", '<svg aria-hidden="true"')
        return svg_content

    # Try to load the icon SVG
    base_path = os.path.join(template_dir, "node_modules", "@material-design-icons", "svg")

    # Try filled version first, then outlined
    for variant in ["filled", "outlined"]:
        icon_path = os.path.join(base_path, variant, f"{icon_name}.svg")
        try:
            with open(icon_path, "r") as f:
                return process_svg(f.read())
        except FileNotFoundError:
            continue

    # If we get here, icon not found in any variant
    lint_error(f"Material Design Icon '{icon_name}' not found")
    return ""


# Export the function so it can be used in templates
template_functions = {"material_icon": get_material_icon}
