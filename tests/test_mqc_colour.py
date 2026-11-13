import pytest

from multiqc.plots.table_object import SECTION_COLORS
from multiqc.utils.mqc_colour import color_to_rgb_string


@pytest.mark.parametrize(
    "color, expected",
    [
        (None, "rgb(0,0,0)"),
        ("rgb(255,0,0)", "rgb(255,0,0)"),
        ("rgba(255,0,0,0.5)", "rgba(255,0,0,0.5)"),
        ("77,175,74", "rgb(77,175,74)"),
        (" 77, 175 ,74 ", "rgb(77,175,74)"),
        ("#4daf4a", "rgb(77,175,74)"),
        ("#ffffff", "rgb(255,255,255)"),
        ("#000000", "rgb(0,0,0)"),
        ("red", "rgb(255,0,0)"),
    ],
)
def test_color_to_rgb_string(color, expected):
    assert color_to_rgb_string(color) == expected


@pytest.mark.parametrize("color", SECTION_COLORS)
def test_table_section_colors_convert_without_warning(color, caplog):
    """Table columns store their color as a bare "r,g,b" triplet, and violin
    plots pass it straight to color_to_rgb_string. Those used to fall back to
    black with an "Error converting color" warning."""
    assert color_to_rgb_string(color) == f"rgb({color})"
    assert "Error converting color" not in caplog.text
