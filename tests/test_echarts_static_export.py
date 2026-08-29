"""
SSR static export tests for the ECharts backend (`multiqc/plots/echarts/static_export.py`).

Uses the same bar plot as `tests/test_plots.py::test_barplot` and
`tests/test_plots_echarts.py`. The autouse `reset` fixture in `tests/conftest.py` resets
`config`/`report` after every test. Skipped entirely when mini-racer/resvg-py are not
installed (they are optional, `pip install 'multiqc[echarts]'`).
"""

import base64
from typing import Dict

import pytest

pytest.importorskip("py_mini_racer")
pytest.importorskip("resvg_py")

from multiqc import config
from multiqc.plots import bargraph, echarts, table, violin
from multiqc.plots.echarts import static_export
from multiqc.plots.table_object import ColumnDict
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


def _make_violin_plot():
    data = {f"Sample{i}": {"metric_a": float(i), "metric_b": float(i * 2), "metric_c": float(i % 5)} for i in range(10)}
    headers: Dict[str, ColumnDict] = {
        "metric_a": {"title": "Metric A"},
        "metric_b": {"title": "Metric B"},
        "metric_c": {"title": "Metric C"},
    }
    return violin.plot(data, headers, table.TableConfig(id="violinplot", title="Test: Violin"))


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


def test_violin_ssr_renders_kde_polygon():
    """
    Proves the `__FN__` renderItem sentinel round-trips through mini-racer: the KDE
    polygon custom series is drawn as an SVG `<polygon>` element, exactly as in
    `scripts/06_violin_final.py` (the reference prototype this builder ports).
    """
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    svg = static_export.render_svg(option, 900, 500)

    assert "<polygon" in svg
