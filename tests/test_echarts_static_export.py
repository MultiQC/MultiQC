"""
SSR static export tests for the ECharts backend (`multiqc/plots/echarts/static_export.py`).

Uses the same bar plot as `tests/test_plots.py::test_barplot` and
`tests/test_plots_echarts.py`. The autouse `reset` fixture in `tests/conftest.py` resets
`config`/`report` after every test. Skipped entirely when mini-racer/resvg-py are not
installed (they are optional, `pip install 'multiqc[echarts]'`).
"""

import base64

import pytest

pytest.importorskip("py_mini_racer")
pytest.importorskip("resvg_py")

from multiqc import config
from multiqc.plots import bargraph, echarts
from multiqc.plots.echarts import static_export
from multiqc.types import Anchor

_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _make_bar_plot():
    return bargraph.plot(
        {
            "Sample0": {},
            "Sample1": {"Cat1": 1},
            "Sample2": {"Cat1": 1, "Cat2": 1},
            "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
        },
        ["Cat1", "Cat2"],
        bargraph.BarPlotConfig(id="bargraph", title="Test: Bar Graph"),
    )


def test_flat_plot_embeds_base64_png():
    config.plotting_engine = "echarts"
    config.plots_force_flat = True

    plot = _make_bar_plot()
    html = plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    marker = '<img src="data:image/png;base64,'
    assert marker in html

    start = html.index(marker) + len(marker)
    end = html.index('"', start)
    png_bytes = base64.b64decode(html[start:end])

    assert png_bytes[:8] == _PNG_MAGIC
    assert len(png_bytes) > 0


def test_png_rendering_is_deterministic():
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    svg1 = static_export.render_svg(option, 1100, 500)
    svg2 = static_export.render_svg(option, 1100, 500)

    png1 = static_export.svg_to_png(svg1)
    png2 = static_export.svg_to_png(svg2)

    assert png1 == png2
    assert png1[:8] == _PNG_MAGIC
    assert len(png1) > 0
