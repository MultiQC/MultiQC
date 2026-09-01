"""
Server-side static image export for the ECharts plotting engine.

Renders ECharts options to SVG in-process, via a MiniRacer (V8) context running the
SAME committed browser bundle used by the interactive template, then rasterises to PNG
via resvg. No subprocess/batch machinery: mini-racer is in-process and fast (tens of
milliseconds per plot), unlike the batch-export workarounds the old Plotly static path
needed to work around its export engine's instability.

See "Static export path (SSR)" in `multiqc-echarts-exploration/BUILD_PLAN.md` for the
contract this module implements.
"""

import base64
import importlib.resources
import json
import logging
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from multiqc import config
from multiqc.core import tmp_dir
from multiqc.plots.echarts import theme
from multiqc.plots.layout import LayoutIR
from multiqc.utils.util_functions import dump_json

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

logger = logging.getLogger(__name__)

_BUNDLE_PACKAGE = "multiqc.templates.default"
_BUNDLE_RESOURCE = "assets/js/packages/echarts-6.1.0.custom.min.js"
_FONT_PACKAGE = "multiqc.plots.echarts"
_FONT_RESOURCE = "fonts/DejaVuSans.ttf"

_INSTALL_MSG = "Static ECharts export requires mini-racer and resvg-py: pip install 'multiqc[echarts]'"

# `__FN__` sentinels (`{"__FN__": True, "body": "<js source>"}`) are the only way a real JS
# function crosses the Python/JS bridge (only the violin builder emits them, Phase 2). The
# bodies are MultiQC-generated string constants, never user data, so `new Function(...)` on
# them is safe; this walker generalizes the pattern from `scripts/06_violin_final.py:20-43`
# to every series' `renderItem` and every axis' `axisLabel.formatter`.
_FN_WALKER_JS = """
function __mqcToFn(sentinel, argNames) {
    return Function.apply(null, argNames.concat([sentinel.body]));
}
function __mqcApplySentinels(opt) {
    (opt.series || []).forEach(function (s) {
        if (s.renderItem && s.renderItem.__FN__) {
            s.renderItem = __mqcToFn(s.renderItem, ["params", "api"]);
        }
    });
    ["xAxis", "yAxis"].forEach(function (axisKey) {
        var axes = opt[axisKey];
        if (!axes) return;
        (Array.isArray(axes) ? axes : [axes]).forEach(function (axis) {
            var formatter = axis.axisLabel && axis.axisLabel.formatter;
            if (formatter && formatter.__FN__) {
                axis.axisLabel.formatter = __mqcToFn(formatter, ["v"]);
            }
        });
    });
}
"""

_RENDER_JS = """
function render(optionJson, w, h) {
    var opt = JSON.parse(optionJson);
    __mqcApplySentinels(opt);
    var chart = echarts.init(null, "multiqc-light", { renderer: "svg", ssr: true, width: w, height: h });
    chart.setOption(opt);
    var svg = chart.renderToSVGString();
    chart.dispose();
    return svg;
}
"""

_racer_instance: Optional[Any] = None

_PDF_INSTALL_MSG = "PDF export requires the pdf extra: pip install 'multiqc[pdf]'"


def _racer() -> Any:
    """
    Lazily create the module-level MiniRacer singleton: one V8 context per process,
    bootstrapped once with the vendored ECharts bundle, the light SSR theme, and the
    `render(optionJson, w, h)` entry point.
    """
    global _racer_instance
    if _racer_instance is not None:
        return _racer_instance

    try:
        from py_mini_racer import MiniRacer
    except ImportError as e:
        raise RuntimeError(_INSTALL_MSG) from e

    bundle_js = importlib.resources.files(_BUNDLE_PACKAGE).joinpath(_BUNDLE_RESOURCE).read_text(encoding="utf-8")

    bootstrap = "\n".join(
        [
            # Alias the globals the UMD/IIFE bundle expects when there is no browser window.
            "var self = this; var window = this; var global = this;",
            bundle_js,
            f'echarts.registerTheme("multiqc-light", {json.dumps(theme.LIGHT_THEME)});',
            _FN_WALKER_JS,
            _RENDER_JS,
        ]
    )

    ctx = MiniRacer()
    ctx.eval(bootstrap)
    _racer_instance = ctx
    return ctx


def render_svg(option: Dict[str, Any], width: int, height: int) -> str:
    """
    Render one ECharts option to an SVG string via the in-process SSR bundle.

    Uses `dump_json` (not plain `json.dumps`): dataset values can legitimately contain
    NaN (e.g. a bar category missing for one sample), and Python's `json.dumps` emits the
    non-standard `NaN` token, which a real JSON parser (the bundle's `JSON.parse`) rejects.
    `dump_json` recursively replaces NaN/Infinity with `null`, same as the interactive
    payload embedded in the HTML report.
    """
    ctx = _racer()
    svg = ctx.call("render", dump_json(option), width, height)
    return str(svg)


def svg_to_png(svg: str) -> bytes:
    """
    Rasterise an SVG string to PNG bytes via resvg. Font is pinned to a single family
    (`theme.FONT_FAMILY`) and shipped as a ttf so rendering is identical on every machine
    (`skip_system_fonts=True` disables substitution from whatever happens to be installed).
    `zoom=2` mirrors the Plotly flat-export retina scale (`plot.py`'s `scale=2.0`).
    """
    try:
        import resvg_py
    except ImportError as e:
        raise RuntimeError(_INSTALL_MSG) from e

    font_traversable = importlib.resources.files(_FONT_PACKAGE).joinpath(_FONT_RESOURCE)
    with importlib.resources.as_file(font_traversable) as font_path:
        png = resvg_py.svg_to_bytes(
            svg_string=svg,
            font_files=[str(font_path)],
            skip_system_fonts=True,
            zoom=2,
        )
    return bytes(png)


def _svg_to_pdf(svg: str) -> bytes:
    """
    Convert a rendered SVG to a VECTOR PDF via svglib + reportlab (both pure-Python, no
    system libraries). Lazily imported: the converters ship only as the optional `pdf`
    extra, so PDF requested without them installed raises a clear, actionable error rather
    than failing obscurely.

    One known limitation: reportlab renders a many-stop linear gradient with only its end
    stops, so a heatmap's `visualMap` colour-bar legend degrades to a two-tone bar. The
    heatmap cells and their printed values are unaffected (they carry solid per-cell
    colours), and every other plot type is faithful. SVG export (also vector) keeps the
    full gradient if that legend matters.
    """
    try:
        import io

        from reportlab.graphics import renderPDF  # type: ignore[import-untyped]
        from svglib.svglib import svg2rlg  # type: ignore[import-untyped]
    except ImportError as e:
        raise RuntimeError(_PDF_INSTALL_MSG) from e

    drawing = svg2rlg(io.BytesIO(svg.encode("utf-8")))
    if drawing is None:
        raise RuntimeError("svglib could not parse the rendered plot SVG for PDF export")
    buf = io.BytesIO()
    renderPDF.drawToFile(drawing, buf, showBoundary=0)
    return buf.getvalue()


def _dataset_height(plot: "Plot[Any, Any]", dataset: Any) -> int:
    """
    The pixel height for one dataset's rendered chart: `plot.layout_ir` merged with the
    per-dataset `dataset.layout` fragment, same merge `Plot.get_figure` performs for the
    Plotly path. Width is not read from here: flat plots force it to 600/1100 (see
    `flat_plot_html`), matching `Plot.get_figure`'s `flat=True` branch.
    """
    merged = plot.layout_ir.merged_with(LayoutIR.from_dataset_layout(dataset.layout))
    if merged.height is None:
        raise ValueError(f"Plot {plot.id!r} has no layout height set")
    return int(merged.height)


def _write_export_formats(svg: str, file_name: str) -> Optional[bytes]:
    """
    When `config.export_plots`, write every format in `config.export_plot_formats` under
    `tmp_dir.plots_tmp_dir()`. Returns the rendered PNG bytes if a PNG was produced (so
    callers embedding a PNG in the HTML can reuse it instead of re-rendering).
    """
    png_bytes: Optional[bytes] = None
    for file_ext in config.export_plot_formats:
        if file_ext == "pdf":
            content: bytes = _svg_to_pdf(svg)
        elif file_ext == "svg":
            content = svg.encode("utf-8")
        elif file_ext == "png":
            png_bytes = svg_to_png(svg)
            content = png_bytes
        else:
            continue

        plot_path = tmp_dir.plots_tmp_dir() / file_ext / f"{file_name}.{file_ext}"
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        plot_path.write_bytes(content)

    return png_bytes


def _render_variant(
    plot: "Plot[Any, Any]",
    ds_idx: int,
    is_pct: bool,
    is_log: bool,
    active: bool,
    file_name: str,
    width: int,
    embed_in_html: bool,
    plots_dir_name: Optional[str],
) -> str:
    from multiqc.plots import echarts as echarts_pkg  # lazy: avoids a circular import

    dataset = plot.datasets[ds_idx]
    height = _dataset_height(plot, dataset)
    option = echarts_pkg.get_option(plot, ds_idx, is_log=is_log, is_pct=is_pct)
    svg = render_svg(option, width, height)

    png_bytes: Optional[bytes] = None
    if config.export_plots:
        png_bytes = _write_export_formats(svg, file_name)

    if embed_in_html:
        if png_bytes is None:
            png_bytes = svg_to_png(svg)
        img_src = f"data:image/png;base64,{base64.b64encode(png_bytes).decode('utf8')}"
    else:
        if plots_dir_name is None:
            raise ValueError("plots_dir_name is required for non-embedded plots")
        if png_bytes is None:
            # PNG is always needed for the <img> tag even when config.export_plots is off
            # or "png" isn't in config.export_plot_formats.
            png_bytes = svg_to_png(svg)
            png_path = tmp_dir.plots_tmp_dir() / "png" / f"{file_name}.png"
            png_path.parent.mkdir(parents=True, exist_ok=True)
            png_path.write_bytes(png_bytes)
        img_src = str(Path(plots_dir_name) / "png" / f"{file_name}.png")

    style = "" if active else "display:none;"
    return "".join(
        [
            f'<div class="mqc_mplplot" style="{style}" id="{file_name}">',
            f'<img src="{img_src}" height="{height}px" width="{width}px"/>',
            "</div>",
        ]
    )


def flat_plot_html(plot: "Plot[Any, Any]", embed_in_html: bool, plots_dir_name: Optional[str]) -> str:
    """
    Mirrors the `figures_to_export` loop of `Plot.flat_plot`: per dataset, the same
    `-cnt` / `-pct` / `-log` / `-pct-log` variants, with identical div ids and active-flag
    logic (the flat-toggle JS in `src/js/flat.js` keys on those ids), each rendered via SSR
    instead of the old Plotly static exporter.

    Every plot type is supported; an unknown one raises `NotImplementedError` (via
    `echarts.get_option` in `_render_variant`) rather than rendering a placeholder.
    """
    if not embed_in_html and plots_dir_name is None:
        raise ValueError("plots_dir_name is required for non-embedded plots")

    width = 600 if config.simple_output else 1100

    # (ds_idx, is_pct, is_log, active, file_name), same variant logic as Plot.flat_plot.
    variants: List[Any] = []
    for ds_idx, dataset in enumerate(plot.datasets):
        variants.append(
            (
                ds_idx,
                False,
                False,
                ds_idx == 0 and not plot.p_active and not plot.l_active,
                dataset.uid if not plot.add_log_tab and not plot.add_pct_tab else f"{dataset.uid}-cnt",
            )
        )
        if plot.add_pct_tab:
            variants.append((ds_idx, True, False, ds_idx == 0 and plot.p_active, f"{dataset.uid}-pct"))
        if plot.add_log_tab:
            variants.append((ds_idx, False, True, ds_idx == 0 and plot.l_active, f"{dataset.uid}-log"))
        if plot.add_pct_tab and plot.add_log_tab:
            variants.append(
                (ds_idx, True, True, ds_idx == 0 and plot.p_active and plot.l_active, f"{dataset.uid}-pct-log")
            )

    return "".join(
        _render_variant(
            plot,
            ds_idx,
            is_pct,
            is_log,
            active,
            file_name,
            width,
            embed_in_html,
            plots_dir_name,
        )
        for ds_idx, is_pct, is_log, active, file_name in variants
    )


def interactive_html(plot: "Plot[Any, Any]", ds_idx: int) -> str:
    """
    Build a self-contained interactive ECharts `<div>` + `<script>` HTML fragment for one
    dataset: the same option (`echarts.get_option`) and browser bundle the report's own
    ECharts template renders, so `Plot.show()`/`Plot.save()` in a notebook match the
    report when `config.plotting_engine == "echarts"`.

    Unlike the report, this fragment carries the ECharts bundle and the `__FN__` sentinel
    revival inline, so it works standalone in a notebook cell or a bare saved HTML file
    with no surrounding report scaffolding. The bundle is guarded behind
    `typeof echarts === "undefined"` so several plots rendered in one notebook do not each
    re-register the theme/bundle redundantly, though each still inlines its own copy of the
    script (unavoidable without a shared asset server).
    """
    from multiqc.plots import echarts as echarts_pkg  # lazy: avoids a circular import

    dataset = plot.datasets[ds_idx]
    width = 600 if config.simple_output else 1100
    height = _dataset_height(plot, dataset)
    option = echarts_pkg.get_option(plot, ds_idx, is_log=plot.l_active, is_pct=plot.p_active)

    bundle_js = importlib.resources.files(_BUNDLE_PACKAGE).joinpath(_BUNDLE_RESOURCE).read_text(encoding="utf-8")
    div_id = f"mqc-nb-{plot.id}-{uuid.uuid4().hex[:8]}"

    script = "\n".join(
        [
            f'<script id="{div_id}-script">',
            'if (typeof echarts === "undefined") {',
            bundle_js,
            "}",
            f'echarts.registerTheme("multiqc-light", {json.dumps(theme.LIGHT_THEME)});',
            _FN_WALKER_JS,
            f"var __mqcOption = JSON.parse({json.dumps(dump_json(option))});",
            "__mqcApplySentinels(__mqcOption);",
            f'echarts.init(document.getElementById("{div_id}"), "multiqc-light", '
            f"{{width: {width}, height: {height}}}).setOption(__mqcOption);",
            "</script>",
        ]
    )
    return f'<div id="{div_id}" style="width:{width}px;height:{height}px;"></div>\n{script}'
