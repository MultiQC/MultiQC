"""
Serialization contract tests for the ECharts backend (`multiqc/plots/echarts/`).

Uses the same bar plot as `tests/test_plots.py::test_barplot`. The autouse `reset`
fixture in `tests/conftest.py` resets `config`/`report` after every test.
"""

import json
import math
from typing import Any, Dict, List, Union

import multiqc
from multiqc import config, report
from multiqc.core.update_config import ClConfig
from multiqc.plots import bargraph, box, echarts, heatmap, linegraph, scatter, seqcontent, table, violin
from multiqc.plots.table_object import ColumnDict
from multiqc.types import Anchor, PlotType


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


def _make_grouped_bar_plot():
    # Mirrors the real-world usage in multiqc/modules/ribowaltz/ribowaltz.py: one
    # sample_groups group per region, each containing one row per sample; rows for the
    # same sample share an offset_group id ("sample1"/"sample2") so they align across
    # regions.
    data = {
        "sample1_5utr": {"Frame 0": 10, "Frame 1": 5, "Frame 2": 2},
        "sample2_5utr": {"Frame 0": 8, "Frame 1": 4, "Frame 2": 1},
        "sample1_cds": {"Frame 0": 20, "Frame 1": 3, "Frame 2": 1},
        "sample2_cds": {"Frame 0": 18, "Frame 1": 2, "Frame 2": 1},
    }
    sample_groups = {
        "5utr": [["sample1_5utr", "sample1"], ["sample2_5utr", "sample2"]],
        "cds": [["sample1_cds", "sample1"], ["sample2_cds", "sample2"]],
    }
    return bargraph.plot(
        data,
        ["Frame 0", "Frame 1", "Frame 2"],
        bargraph.BarPlotConfig(
            id="bargraph_grouped",
            title="Test: Grouped Bar Graph",
            sample_groups=sample_groups,
        ),
    )


def _make_line_plot():
    return linegraph.plot(
        {
            "Sample1": {0: 1.0, 1: 2.0},
            "Sample2": {0: 3.0, 1: 4.0, 2: 5.0},
        },
        linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
    )


def _make_categorical_line_plot():
    return linegraph.plot(
        {
            "Sample1": {"a": 1.0, "b": 2.0},
            "Sample2": {"a": 3.0, "b": 4.0},
        },
        linegraph.LinePlotConfig(id="linegraph_categorical", title="Test: Categorical Line Graph", categories=True),
    )


def _make_line_plot_with_bands():
    return linegraph.plot(
        {"Sample1": {0: 1.0, 1: 2.0}},
        linegraph.LinePlotConfig(
            id="linegraph_bands",
            title="Test: Line Graph With Bands",
            y_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
        ),
    )


def _make_scatter_plot():
    return scatter.plot(
        {
            "Sample1": [{"x": 1, "y": 2, "color": "rgb(255,0,0)", "marker_size": 12}],
            "Sample2": [{"x": 3, "y": 4}],
        },
        scatter.ScatterConfig(id="scatterplot", title="Test: Scatter"),
    )


def _make_bar_plot_with_bands():
    return bargraph.plot(
        {"Sample1": {"Cat1": 1}, "Sample2": {"Cat1": 2}},
        ["Cat1"],
        bargraph.BarPlotConfig(
            id="bargraph_bands",
            title="Test: Bar Graph With Bands",
            x_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
            x_lines=[{"value": 1.5, "label": "threshold"}],
            # Meaningless for bar's category (sample) yAxis: must be dropped, not crash.
            y_bands=[{"from": 0, "to": 1, "color": "#ff0000"}],
            y_lines=[{"value": 0.5, "label": "y_ignored"}],
        ),
    )


def _make_scatter_plot_with_bands():
    return scatter.plot(
        {"Sample1": [{"x": 1, "y": 2}], "Sample2": [{"x": 3, "y": 4}]},
        scatter.ScatterConfig(
            id="scatterplot_bands",
            title="Test: Scatter With Bands",
            x_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
            y_lines=[{"value": 2.5, "label": "y_threshold"}],
        ),
    )


def _make_large_scatter_plot(n_points: int):
    return scatter.plot(
        {f"Sample{i}": [{"x": i, "y": i}] for i in range(n_points)},
        scatter.ScatterConfig(id="scatterplot_large", title="Test: Large Scatter"),
    )


def _make_heatmap_plot(display_values=False):
    return heatmap.plot(
        data=[[1, 2], [3, 4]],
        xcats=["Cat1", "Cat2"],
        ycats=["Sample1", "Sample2"],
        pconfig=heatmap.HeatmapConfig(id="echarts_heatmap", title="Test: Heatmap", display_values=display_values),
    )


def _make_clustered_heatmap_plot():
    # Same fixture as test_plots.py::test_heatmap_clustering_produces_reordered_data:
    # designed so clustering actually reorders rows/cats (a 2x2 or uniform matrix may not).
    data = [
        [10, 0, 11],
        [0, 10, 0],
        [11, 0, 10],
    ]
    return heatmap.plot(
        data=data,
        xcats=["X1", "X2", "X3"],
        ycats=["A", "B", "C"],
        pconfig=heatmap.HeatmapConfig(
            id="echarts_heatmap_clustered",
            title="Test: Clustered Heatmap",
            cluster_rows=True,
            cluster_cols=True,
        ),
    )


def _make_box_plot():
    return box.plot(
        {"Sample1": [1.0, 2.0, 3.0, 4.0], "Sample2": [5.0, 6.0, 7.0, 8.0]},
        box.BoxPlotConfig(id="boxplot", title="Test: Box"),
    )


def _make_stats_box_plot():
    return box.plot(
        {
            "Sample1": {"min": 1.0, "q1": 2.0, "median": 3.0, "q3": 4.0, "max": 5.0},
            "Sample2": {"min": 2.0, "q1": 3.0, "median": 4.0, "q3": 5.0, "max": 6.0},
        },
        box.BoxPlotConfig(id="boxplot_stats", title="Test: Box Stats"),
    )


def _make_violin_plot():
    data = {f"Sample{i}": {"metric_a": float(i), "metric_b": float(i * 2), "metric_c": float(i % 5)} for i in range(10)}
    headers: Dict[str, ColumnDict] = {
        "metric_a": {"title": "Metric A"},
        "metric_b": {"title": "Metric B"},
        "metric_c": {"title": "Metric C"},
    }
    return violin.plot(data, headers, table.TableConfig(id="violinplot", title="Test: Violin"))


def test_interactive_plot_adds_echarts_key_when_engine_is_echarts():
    config.plotting_engine = "echarts"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["yAxis"]["type"] == "category"
    assert "series" not in skeleton


def test_interactive_plot_title_strips_html():
    config.plotting_engine = "echarts"
    plot = bargraph.plot(
        {
            "Sample1": {"Cat1": 1},
        },
        ["Cat1"],
        bargraph.BarPlotConfig(id="bargraph_title", title="Foo Bar<br><sup>42 things</sup>"),
    )
    assert isinstance(plot, bargraph.BarPlot)  # narrow BarPlot | str | None for mypy
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["title"]["text"] == "Foo Bar"
    assert skeleton["title"]["subtext"] == "42 things"


def test_interactive_plot_unsupported_plot_type_does_not_raise(monkeypatch):
    # All PlotTypes with a figure are ported now (bar/line/scatter/heatmap/box/violin);
    # exercise the generic "unsupported" fallback path directly instead, by pretending
    # box plots aren't ported, so a single unported type still can't crash the report.
    monkeypatch.delitem(echarts._BUILDERS, PlotType.BOX)
    config.plotting_engine = "echarts"
    plot = _make_box_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert dumped["echarts"] == {"unsupported": "box plot"}


def test_interactive_plot_adds_echarts_key_for_line_plot():
    config.plotting_engine = "echarts"
    plot = _make_line_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "value"
    # No slider (matches Plotly, which has no mini-plot strip); Plotly-style click+drag
    # box-zoom on both axes instead (POLISH.md #17), via the toolbox dataZoom feature
    # plus one "inside" dataZoom per zoomable axis holding the current zoom range, with
    # mouse-wheel zoom disabled so it doesn't hijack page scrolling. See
    # `multiqc/plots/echarts/converter.py::zoom_option` and `line.py::layout_option`.
    assert skeleton["dataZoom"] == [
        {"type": "inside", "xAxisIndex": [0], "zoomOnMouseWheel": False, "moveOnMouseWheel": False},
        {"type": "inside", "yAxisIndex": [0], "zoomOnMouseWheel": False, "moveOnMouseWheel": False},
    ]
    assert skeleton["toolbox"]["feature"]["dataZoom"]["xAxisIndex"] == [0]
    assert skeleton["toolbox"]["feature"]["dataZoom"]["yAxisIndex"] == [0]
    assert "series" not in skeleton


def test_interactive_plot_includes_echarts_key_by_default():
    # ECharts is the default engine now, so the serialized payload carries the echarts key.
    assert config.plotting_engine == "echarts"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))
    assert "echarts" in report.plot_data[plot.anchor]


def test_interactive_plot_omits_echarts_key_for_plotly_engine():
    config.plotting_engine = "plotly"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))
    assert "echarts" not in report.plot_data[plot.anchor]


def test_get_option_series_length_matches_categories():
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == len(dataset.cats)
    for cat, series_option in zip(dataset.cats, option["series"]):
        assert series_option["data"] == cat.data


def test_get_option_is_pct_uses_data_pct():
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=True)

    dataset = plot.datasets[0]
    for cat, series_option in zip(dataset.cats, option["series"]):
        # NaN != NaN, so compare positionally while tolerating NaN entries.
        for actual, expected in zip(series_option["data"], cat.data_pct):
            if math.isnan(expected):
                assert math.isnan(actual)
            else:
                assert actual == expected


def test_get_option_grouped_bar_produces_stacked_series_per_sample():
    # FIX-NEEDED #3 in multiqc-echarts-exploration/PARITY.md: grouped bars
    # (`dataset.group_labels` set, e.g. by `pconfig.sample_groups`) used to throw.
    plot = _make_grouped_bar_plot()
    dataset = plot.datasets[0]
    assert dataset.group_labels is not None

    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    unique_groups = list(dict.fromkeys(dataset.group_labels))
    assert len(unique_groups) == 2  # "5utr" and "cds"
    assert option["yAxis"]["data"] == unique_groups

    # One series per (unique sample, category): ECharts has no native multicategory
    # axis, so each sample becomes a `stack` id positioned at its group's slot instead
    # of Plotly's shared `offsetgroup`.
    unique_samples = list(dict.fromkeys(dataset.samples))
    assert len(option["series"]) == len(unique_samples) * len(dataset.cats)

    stack_ids = {s["stack"] for s in option["series"]}
    assert dataset.offset_groups is not None
    assert stack_ids == set(dataset.offset_groups.values())

    for s in option["series"]:
        assert len(s["data"]) == len(unique_groups)
        # Every row has a value in exactly one group slot for this dataset (each
        # sample appears once per region); the rest are None placeholders.
        assert sum(1 for v in s["data"] if v is not None) == 1


def test_get_option_is_log_sets_switch_controlled_axis_type():
    plot = _make_bar_plot()
    assert plot.axis_controlled_by_switches == ["xaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["xAxis"]["type"] == "log"


def test_get_option_bar_value_axis_is_not_scaled():
    # Bars must anchor at 0 (matches Plotly, which pins an explicit `minallowed=0` for
    # bar plots, see `bargraph.py`); `scale: true` is for line/scatter/box only.
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["min"] == 0
    assert "scale" not in option["xAxis"]


def test_get_option_value_axis_has_si_formatter_sentinel_category_axis_does_not():
    # POLISH.md #12: value axes get an SI-abbreviation `__FN__` sentinel (same pattern
    # violin.py's yAxis uses); bar's yAxis is a category axis (sample names) and must
    # stay untouched, matching Plotly (no formatter needed on category ticks).
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["type"] == "value"
    assert option["xAxis"]["axisLabel"]["formatter"]["__FN__"] is True
    assert option["yAxis"]["type"] == "category"
    assert "formatter" not in option["yAxis"].get("axisLabel", {})


def test_si_axis_formatter_body_golden_values():
    """
    GOLDEN cross-language contract (POLISH.md #12): `_si_axis_formatter_body`'s JS source
    must stay in lockstep with `formatAxisNumber()` in
    `multiqc/templates/default/src/js/echarts-plotting.js` (same duplication pattern as
    violin.py's kde()/JS kde() pair). Executed here in a real JS engine (MiniRacer,
    already a project dependency via `static_export.py`) to assert it SI-abbreviates the
    way Plotly's own axis ticks do: ~3 significant figures, trailing zeros trimmed.
    """
    from py_mini_racer import MiniRacer

    from multiqc.plots.echarts.converter import _si_axis_formatter_body

    ctx = MiniRacer()
    ctx.eval(f"function fmt(v) {{ {_si_axis_formatter_body('')} }}")
    assert ctx.call("fmt", 450000000) == "450M"
    assert ctx.call("fmt", 12000) == "12k"
    assert ctx.call("fmt", 1500000) == "1.5M"
    assert ctx.call("fmt", 2000000000) == "2G"
    assert ctx.call("fmt", 20) == "20"
    assert ctx.call("fmt", 0) == "0"
    assert ctx.call("fmt", -5000) == "-5k"

    # Suffix (e.g. a coverage "x" unit) combines with the SI abbreviation instead of the
    # two fighting over the same tick.
    ctx.eval(f"function fmtx(v) {{ {_si_axis_formatter_body('x')} }}")
    assert ctx.call("fmtx", 12000) == "12kx"


def test_get_option_line_series_matches_lines():
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == len(dataset.lines)
    for line, series_option in zip(dataset.lines, option["series"]):
        assert series_option["type"] == "line"
        assert series_option["name"] == line.name
        assert series_option["data"] == [list(p) for p in line.pairs]
        assert series_option["smooth"] is False


def test_get_option_line_value_axis_has_no_static_data():
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "data" not in option["xAxis"]


def test_get_option_line_categorical_axis_uses_plain_values():
    plot = _make_categorical_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["type"] == "category"
    assert option["xAxis"]["data"] == ["a", "b"]

    dataset = plot.datasets[0]
    for line, series_option in zip(dataset.lines, option["series"]):
        assert series_option["data"] == [y for _, y in line.pairs]


def test_get_option_line_value_axes_are_scaled():
    # Non-categorical x-axis and y-axis are both value axes for a line plot: `scale:
    # true` fits them to the data extent instead of ECharts' default of forcing 0 in,
    # matching Plotly's autorange (POLISH.md #5).
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["scale"] is True


def test_get_option_line_categorical_axis_is_not_scaled():
    # A categorical x-axis has no meaningful "scale to data" behavior; only the
    # (still value-typed) y-axis gets it.
    plot = _make_categorical_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "scale" not in option["xAxis"]
    assert option["yAxis"]["scale"] is True


def test_get_option_line_is_log_sets_yaxis():
    plot = _make_line_plot()
    assert plot.axis_controlled_by_switches == ["yaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["yAxis"]["type"] == "log"


def test_get_option_line_y_bands_produce_silent_markarea_series():
    plot = _make_line_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    assert marker_series[0]["markArea"]["data"] == [
        [{"yAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"yAxis": 1}]
    ]

    # The skeleton (interactive JS path) carries the same payload under `_mqc.bandsLines`.
    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markArea"] == marker_series[0]["markArea"]


def test_get_option_bar_x_bands_and_lines_produce_silent_marker_series():
    plot = _make_bar_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    # x_bands/x_lines land on the value axis (ECharts "xAxis"): bar is horizontal, so
    # xAxis is the value axis (see converter.bands_and_lines docstring).
    assert marker_series[0]["markArea"]["data"] == [
        [{"xAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"xAxis": 1}]
    ]
    assert marker_series[0]["markLine"]["data"][0]["xAxis"] == 1.5

    # y_bands/y_lines target bar's category (sample) axis, which is meaningless for a
    # numeric threshold: dropped cleanly, not crashed and not emitted.
    assert "yAxis" not in marker_series[0]["markArea"]["data"][0][0]
    assert all("yAxis" not in entry for entry in marker_series[0]["markLine"]["data"])

    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markArea"] == marker_series[0]["markArea"]


def test_get_option_scatter_x_bands_and_y_lines_produce_silent_marker_series():
    plot = _make_scatter_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    assert marker_series[0]["markArea"]["data"] == [
        [{"xAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"xAxis": 1}]
    ]
    assert marker_series[0]["markLine"]["data"] == [
        {
            "yAxis": 2.5,
            "lineStyle": {"color": None, "width": 2, "type": "solid"},
            "label": {"formatter": "y_threshold"},
        }
    ]

    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markLine"] == marker_series[0]["markLine"]


def test_interactive_plot_adds_echarts_key_for_scatter_plot():
    config.plotting_engine = "echarts"
    plot = _make_scatter_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["tooltip"]["trigger"] == "item"
    assert "series" not in skeleton


def test_get_option_scatter_value_axes_are_scaled():
    # Scatter has no category axis: both x and y are value axes, so both get `scale:
    # true` (Plotly-style data-fitted autorange, not forced-0), matching line/box.
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["scale"] is True


def test_get_option_scatter_series_has_one_series_with_value_items():
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == 1
    series_option = option["series"][0]
    assert series_option["type"] == "scatter"
    assert len(series_option["data"]) == len(dataset.points)

    for point, item in zip(dataset.points, series_option["data"]):
        assert item["value"] == [point["x"], point["y"]]
        assert item["name"] == point["name"]

    # Sample1's explicit color/marker_size are carried through per-item.
    sample1_item = series_option["data"][0]
    assert sample1_item["symbolSize"] == 12
    assert sample1_item["itemStyle"]["color"] == "rgb(255,0,0)"


def test_get_option_scatter_axis_has_no_static_data():
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "data" not in option["xAxis"]
    assert "data" not in option["yAxis"]


def test_scatter_mark_count_is_point_count():
    plot = _make_scatter_plot()
    assert echarts.scatter.mark_count(plot.datasets[0]) == 2


def test_serialize_scatter_uses_svg_renderer_below_threshold():
    config.plotting_engine = "echarts"
    config.echarts_canvas_threshold = 3000
    plot = _make_large_scatter_plot(10)
    result = echarts.serialize(plot)
    assert result["renderer"] == "svg"


def test_serialize_scatter_uses_canvas_renderer_above_threshold():
    config.plotting_engine = "echarts"
    config.echarts_canvas_threshold = 3000
    plot = _make_large_scatter_plot(3001)
    result = echarts.serialize(plot)
    assert result["renderer"] == "canvas"


def test_interactive_plot_adds_echarts_key_for_heatmap_plot():
    config.plotting_engine = "echarts"
    plot = _make_heatmap_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "category"
    assert skeleton["yAxis"]["type"] == "category"
    assert "series" not in skeleton
    assert "data" not in skeleton["xAxis"]
    assert "data" not in skeleton["yAxis"]

    # visualMap: converted from the Plotly colorscale stop list, resampled into an evenly
    # spaced colour list ECharts' own even-spacing-only visualMap can consume (see
    # `test_colorscale_colors_resamples_uneven_stops_to_exact_positions` for the resample
    # itself, which is what makes this reproduce Plotly's positions rather than just its
    # ordered colours).
    visual_map = skeleton["visualMap"]
    assert visual_map["calculable"] is True
    assert isinstance(visual_map["inRange"]["color"], list)
    assert len(visual_map["inRange"]["color"]) > 0
    assert visual_map["min"] == 1
    assert visual_map["max"] == 4


def test_get_option_heatmap_builds_xyz_cell_data():
    plot = _make_heatmap_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["data"] == ["Cat1", "Cat2"]
    assert option["yAxis"]["data"] == ["Sample1", "Sample2"]

    series_option = option["series"][0]
    assert series_option["type"] == "heatmap"
    assert sorted(series_option["data"]) == [[0, 0, 1], [0, 1, 3], [1, 0, 2], [1, 1, 4]]


def test_heatmap_mark_count_is_rows_times_cols():
    plot = _make_heatmap_plot()
    assert echarts.heatmap.mark_count(plot.datasets[0]) == 4


def test_get_option_heatmap_cell_labels_when_display_values_enabled():
    plot = _make_heatmap_plot(display_values=True)
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    series_option = option["series"][0]
    for item in series_option["data"]:
        assert isinstance(item, dict)
        assert item["label"]["show"] is True
        _xi, _yi, val = item["value"]
        assert item["label"]["formatter"] == f"{val:.2f}"


def _colorscale_colors_for(colstops, reverse_colors=False):
    plot = heatmap.plot(
        data=[[0]],
        xcats=["Cat1"],
        ycats=["Sample1"],
        pconfig=heatmap.HeatmapConfig(
            id="echarts_heatmap_colorscale",
            title="Test: Heatmap Colorscale",
            colstops=colstops,
            reverse_colors=reverse_colors,
        ),
    )
    assert isinstance(plot, heatmap.HeatmapPlot)  # narrow HeatmapPlot | str | None for mypy
    return echarts.heatmap._colorscale_colors(plot.datasets[0])


def test_colorscale_colors_resamples_uneven_stops_to_exact_positions():
    # FastQC's Status Checks heatmap uses stops that are NOT evenly spaced
    # (0, 0.25, 0.5, 1): naively spreading the raw 4-colour list evenly across [0, 1]
    # (the pre-fix behaviour) put the 0.25/0.5 colours in the wrong place. The resample
    # must reproduce Plotly's own piecewise-linear colorscale lookup at each stop.
    colors = _colorscale_colors_for([[0, "#ffffff00"], [0.25, "#d9534f"], [0.5, "#fee391"], [1, "#5cb85c"]])
    n = len(colors)

    def at(pos):
        return colors[round(pos * (n - 1))]

    assert at(0) == "rgba(255, 255, 255, 0.0)"  # not-run: transparent white
    assert at(0.25) == "rgba(217, 83, 79, 1.0)"  # fail: matches Plotly's #d9534f exactly
    assert at(0.5) == "rgba(254, 227, 145, 1.0)"  # warn: matches Plotly's #fee391 exactly
    assert at(1) == "rgba(92, 184, 92, 1.0)"  # pass: matches Plotly's #5cb85c exactly


def test_colorscale_colors_hard_step_is_not_blended():
    # Duplicate/adjacent stop positions are a deliberate hard colour boundary (e.g. a
    # 2-band scale); the resample must reproduce the same sharp jump, not a blend
    # across the two colours flanking it.
    colors = _colorscale_colors_for([[0, "#ff0000"], [0.5, "#ff0000"], [0.5, "#0000ff"], [1, "#0000ff"]])
    n = len(colors)

    def at(pos):
        return colors[round(pos * (n - 1))]

    assert at(0.49) == "rgba(255, 0, 0, 1.0)"
    assert at(0.51) == "rgba(0, 0, 255, 1.0)"


def test_colorscale_colors_reversescale_mirrors_stop_positions():
    # Plotly's own `reversescale` mirrors stop positions (`1 - pos`), not just the
    # colour order; verified against actual `plotly` rendered pixels for this exact,
    # non-symmetric colorscale (see the module docstring in `echarts/heatmap.py`).
    colors = _colorscale_colors_for(
        [[0, "#ffffff"], [0.25, "#d9534f"], [0.5, "#fee391"], [1, "#5cb85c"]],
        reverse_colors=True,
    )
    n = len(colors)

    def at(pos):
        return colors[round(pos * (n - 1))]

    assert at(0) == "rgba(92, 184, 92, 1.0)"
    assert at(0.5) == "rgba(254, 227, 145, 1.0)"
    assert at(1) == "rgba(255, 255, 255, 1.0)"


def test_get_option_heatmap_clustered_switch_uses_clustered_categories():
    plot = _make_clustered_heatmap_plot()
    dataset = plot.datasets[0]
    assert dataset.rows_clustered is not None  # sanity: clustering actually happened

    plot.pconfig.cluster_switch_clustered_active = True
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["data"] == list(dataset.xcats_clustered)
    assert option["yAxis"]["data"] == list(dataset.ycats_clustered)

    plot.pconfig.cluster_switch_clustered_active = False
    option_unclustered = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option_unclustered["xAxis"]["data"] == list(dataset.xcats)
    assert option_unclustered["yAxis"]["data"] == list(dataset.ycats)


def _make_seqcontent_plot():
    # sample1's first bin (t=100) and sample2's bin are the golden RGB fixtures used in
    # tests/test_seqcontent.py::test_rgb_golden_t_100/test_rgb_golden_even_split:
    # t=100 -> rgb(255, 0, 0); a=c=g=t=25 -> rgb(64, 64, 64). sample1's second bin
    # ("2-3", a range) also exercises a variable-width bin (max_bp == 3).
    data_by_sample = {
        "sample1": {
            "1": {"a": 0.0, "c": 0.0, "g": 0.0, "t": 100.0},
            "2-3": {"a": 25.0, "c": 25.0, "g": 25.0, "t": 25.0},
        },
        "sample2": {
            "1": {"a": 25.0, "c": 25.0, "g": 25.0, "t": 25.0},
        },
    }
    return seqcontent.plot(data_by_sample, {"id": "echarts_seqcontent", "title": "Test: Seq Content"})


def test_interactive_plot_adds_echarts_key_for_seqcontent_plot():
    config.plotting_engine = "echarts"
    plot = _make_seqcontent_plot()
    assert isinstance(plot, seqcontent.SeqContentPlot)  # narrow SeqContentPlot | str | None for mypy
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    # Contract: the skeleton is JSON-safe, no live functions (those only appear via the
    # `__FN__` sentinel in the SSR/get_option path, not the interactive dump).
    json.dumps(dumped["echarts"])

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "value"
    assert skeleton["xAxis"]["min"] == 1
    assert skeleton["xAxis"]["max"] == plot.datasets[0].max_bp + 1
    assert skeleton["yAxis"]["type"] == "category"
    assert skeleton["yAxis"]["inverse"] is True
    assert "series" not in skeleton
    assert "data" not in skeleton["yAxis"]
    assert "visualMap" not in skeleton

    # BUG FIX: a bin only partially inside the zoomed window must still be drawn (its
    # visible portion clipped by the grid), not dropped by dataZoom's default "filter"
    # behavior (see seqcontent.py::layout_option). Both the "inside" dataZoom entries
    # and the toolbox box-select dataZoom feature must carry the override.
    assert skeleton["dataZoom"] == [
        {
            "type": "inside",
            "xAxisIndex": [0],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            "filterMode": "none",
        },
        {
            "type": "inside",
            "yAxisIndex": [0],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            "filterMode": "none",
        },
    ]
    assert skeleton["toolbox"]["feature"]["dataZoom"]["filterMode"] == "none"


def test_get_option_seqcontent_builds_one_custom_series_with_golden_rgb_data():
    plot = _make_seqcontent_plot()
    assert isinstance(plot, seqcontent.SeqContentPlot)  # narrow SeqContentPlot | str | None for mypy
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["yAxis"]["data"] == ["sample1", "sample2"]

    series_list = option["series"]
    assert len(series_list) == 1
    series_option = series_list[0]
    assert series_option["type"] == "custom"

    render_item = series_option["renderItem"]
    assert render_item["__FN__"] is True
    assert isinstance(render_item["body"], str)
    assert "rect" in render_item["body"]

    # [start, end, row_idx, r, g, b, opacity], golden RGB values from bin_rgb (see
    # tests/test_seqcontent.py::test_rgb_golden_t_100/test_rgb_golden_even_split); opacity
    # is always 1 in the SSR/get_option path (no toolbox highlight state server-side).
    assert series_option["data"] == [
        [1, 1, 0, 255, 0, 0, 1],
        [2, 3, 0, 64, 64, 64, 1],
        [1, 1, 1, 64, 64, 64, 1],
    ]


def test_seqcontent_axis_data_returns_yaxis_samples():
    plot = _make_seqcontent_plot()
    assert isinstance(plot, seqcontent.SeqContentPlot)  # narrow SeqContentPlot | str | None for mypy
    dataset = plot.datasets[0]
    assert echarts.seqcontent.axis_data(dataset, plot.pconfig) == [("yAxis", ["sample1", "sample2"])]


def test_seqcontent_mark_count_equals_total_bins():
    plot = _make_seqcontent_plot()
    assert isinstance(plot, seqcontent.SeqContentPlot)  # narrow SeqContentPlot | str | None for mypy
    # 2 bins for sample1 + 1 bin for sample2, not n_samples * max_bp.
    assert echarts.seqcontent.mark_count(plot.datasets[0]) == 3


# GOLDEN quartile test: this fixed input + expected five-number/outlier values is the
# cross-language contract asserted here AND mirrored in a comment block at the top of
# `multiqc/templates/default/src/js/plots-echarts/box.js`. The JS `fiveNumberSummary()`/
# `outliers()` port must produce the same output for the same input (checked manually:
# there is no JS unit-test runner in this repo).
_GOLDEN_BOX_VALUES: List[Union[int, float]] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]
_GOLDEN_FIVE_NUMBER_SUMMARY = [1.0, 3.5, 6.0, 8.5, 10.0]
_GOLDEN_OUTLIERS = [100.0]


def test_box_five_number_summary_golden_values():
    assert echarts.box.five_number_summary(_GOLDEN_BOX_VALUES) == _GOLDEN_FIVE_NUMBER_SUMMARY


def test_box_outliers_golden_values():
    assert echarts.box.outliers(_GOLDEN_BOX_VALUES) == _GOLDEN_OUTLIERS


def test_interactive_plot_adds_echarts_key_for_box_plot():
    config.plotting_engine = "echarts"
    plot = _make_box_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["yAxis"]["type"] == "category"
    assert skeleton["yAxis"]["inverse"] is True
    assert skeleton["tooltip"]["trigger"] == "item"
    assert "series" not in skeleton


def test_get_option_box_series_has_five_number_data():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    assert len(boxplot_series["data"]) == len(dataset.samples)
    for values, entry in zip(dataset.data, boxplot_series["data"]):
        assert entry["value"] == echarts.box.five_number_summary(values)
        # Semi-transparent fill + solid border, so the median (drawn via borderColor)
        # is visible against the box body: POLISH.md #6.
        assert entry["itemStyle"]["color"].startswith("rgba(")
        assert entry["itemStyle"]["color"].endswith(",0.5)")
        assert entry["itemStyle"]["borderColor"].startswith("rgb(")


def test_get_option_box_outliers_become_scatter_series():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    scatter_series = [s for s in option["series"] if s["type"] == "scatter"]
    # Neither sample's 4 evenly-spaced values produce a Tukey outlier.
    total_outliers = sum(len(echarts.box.outliers(values)) for values in dataset.data)
    assert total_outliers == 0
    assert scatter_series == []


def test_get_option_box_axis_has_no_static_sample_data():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["yAxis"]["data"] == list(plot.datasets[0].samples)
    assert "data" not in option["xAxis"]


def test_get_option_box_value_axis_is_scaled():
    # Box plots are horizontal: xAxis is the value axis and gets `scale: true`
    # (Plotly-style data-fitted autorange, POLISH.md #7); yAxis is the (inverted)
    # category axis and is untouched by scale.
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["type"] == "category"
    assert "scale" not in option["yAxis"]


def test_get_option_box_stats_data_uses_precomputed_values_directly():
    plot = _make_stats_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert dataset.is_stats_data is True
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    for stats, entry in zip(dataset.data, boxplot_series["data"]):
        assert entry["value"] == [stats["min"], stats["q1"], stats["median"], stats["q3"], stats["max"]]
    # No raw values to derive outliers from.
    assert not any(s["type"] == "scatter" for s in option["series"])


def test_box_mark_count_is_samples_plus_outliers():
    plot = _make_box_plot()
    dataset = plot.datasets[0]
    expected = len(dataset.samples) + sum(len(echarts.box.outliers(values)) for values in dataset.data)
    assert echarts.box.mark_count(dataset) == expected

    outlier_dataset = box.Dataset(
        **{**dataset.__dict__, "samples": ["S1"], "data": [_GOLDEN_BOX_VALUES], "is_stats_data": False}
    )
    assert echarts.box.mark_count(outlier_dataset) == 1 + len(_GOLDEN_OUTLIERS)


# GOLDEN kde() test: this fixed input + expected output is the cross-language contract
# asserted here AND mirrored in a comment block at the top of `multiqc/plots/echarts/violin.py`
# (and, in Task 2.2, at the top of `templates/default/src/js/plots-echarts/violin.js`). The JS
# `kde()` port must reproduce these same values for the same input.
_GOLDEN_KDE_VALUES = [1.0, 2.0, 3.0, 4.0, 5.0]
_GOLDEN_KDE_XS = [1.0, 3.0, 5.0]
_GOLDEN_KDE_DENSITIES = [0.15916497933387785, 0.1802710624663249, 0.15916497933387785]


def test_kde_golden_values():
    assert echarts.violin.kde(_GOLDEN_KDE_VALUES, _GOLDEN_KDE_XS) == _GOLDEN_KDE_DENSITIES


def _renderitem_body(series_option):
    return series_option["renderItem"]["body"]


def _is_kde_series(series_option):
    return series_option["type"] == "custom" and "polygon" in _renderitem_body(series_option)


def test_interactive_plot_adds_echarts_key_for_violin_plot():
    config.plotting_engine = "echarts"
    plot = _make_violin_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    dataset_entry = dumped["echarts"]["datasets"][0]
    skeleton = dataset_entry["layout"]
    assert skeleton["animation"] is False
    assert "series" not in skeleton

    # PLOTLY-STYLE PER-ROW SUBPLOTS: one grid/xAxis/yAxis per visible metric (3 here),
    # not the single shared grid/axis pair every other type gets.
    assert len(skeleton["grid"]) == 3
    assert len(skeleton["xAxis"]) == 3
    assert len(skeleton["yAxis"]) == 3
    for i, (grid, x_axis, y_axis) in enumerate(zip(skeleton["grid"], skeleton["xAxis"], skeleton["yAxis"])):
        assert grid["top"] and grid["height"]  # percentage strings, see _row_geometry
        assert x_axis["type"] == "value"
        assert x_axis["gridIndex"] == i
        assert x_axis["axisLabel"]["formatter"]["__FN__"] is True  # SI-abbreviation sentinel
        assert y_axis["type"] == "value"
        assert y_axis["gridIndex"] == i
        assert y_axis["min"] == -0.5
        assert y_axis["max"] == 0.5
        # ROW ALIGNMENT fix: the yAxis label is fully hidden (contributes nothing to the
        # grid's left inset) rather than carrying the row title as a tick-label formatter
        # (the old VALUE-AXIS TRICK); the title is instead drawn by a native `title`
        # component entry (see test_get_option_violin_title_array_matches_rows).
        assert y_axis["axisLabel"]["show"] is False
        assert "formatter" not in y_axis["axisLabel"]
        # ROW ALIGNMENT fix: every row's grid shares the same left/right and disables
        # both containLabel and ECharts 6's outerBoundsMode, so no row's rect can be
        # grown/shrunk differently than another's (see violin.py's _row_axes docstring).
        assert grid["containLabel"] is False
        assert grid["outerBoundsMode"] == "none"
        assert grid["left"] == skeleton["grid"][0]["left"]
        assert grid["right"] == skeleton["grid"][0]["right"]

    # Per-row click+drag zoom (POLISH.md #8): the toolbox spans every row's xAxisIndex,
    # and there's one `inside` dataZoom per row, no yAxisIndex anywhere (the y-axis is a
    # fake per-row thickness scale, never worth zooming).
    assert skeleton["toolbox"]["feature"]["dataZoom"]["xAxisIndex"] == [0, 1, 2]
    assert skeleton["toolbox"]["feature"]["dataZoom"]["yAxisIndex"] == []
    assert skeleton["dataZoom"] == [
        {
            "type": "inside",
            "xAxisIndex": [0],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            "filterMode": "none",
        },
        {
            "type": "inside",
            "xAxisIndex": [1],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            "filterMode": "none",
        },
        {
            "type": "inside",
            "xAxisIndex": [2],
            "zoomOnMouseWheel": False,
            "moveOnMouseWheel": False,
            "filterMode": "none",
        },
    ]

    # The skeleton itself must stay plain-JSON-safe (the `__FN__` sentinel is a dict, not
    # an actual function): `serialize()` already asserts this globally via `json.dumps`.

    # VIOLIN EXCEPTION: precomputed KDE polygons ride along per-dataset, keyed by metric.
    violins = dataset_entry["violins"]
    assert set(violins.keys()) == {"metric_a", "metric_b", "metric_c"}
    for payload in violins.values():
        assert len(payload["poly"]) == 2 * echarts.violin.N_BINS
        assert len(payload["range"]) == 2
        assert payload["range"][0] < payload["range"][1]


def test_get_option_violin_series_counts():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    n_metrics = len(dataset.metrics)
    assert n_metrics == 3

    kde_series = [s for s in option["series"] if _is_kde_series(s)]
    scatter_series = [s for s in option["series"] if s["type"] == "scatter"]

    assert len(kde_series) == n_metrics
    assert len(scatter_series) == n_metrics
    # No separate min/max-annotation series any more (removed: each row already has its
    # own real-valued x-axis, so a text label repeating the same range was redundant and
    # read as if the beeswarm dots themselves carried on-canvas labels).
    assert len(option["series"]) == 2 * n_metrics


def test_get_option_violin_scatter_items_shape():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    scatter_by_name = {s["name"]: s for s in option["series"] if s["type"] == "scatter"}
    assert set(scatter_by_name) == set(dataset.metrics)

    for metric, series_option in scatter_by_name.items():
        n_values = len(dataset.violin_value_by_sample_by_metric[metric])
        assert len(series_option["data"]) == n_values
        values = list(dataset.violin_value_by_sample_by_metric[metric].values())
        for item in series_option["data"]:
            assert list(item.keys()) == ["value", "name"]
            # Real value (POLISH.md #6: no more per-row 0..1 normalization), constant
            # y=0 (the row's own vertical center; jitter spreads points visually).
            real_value, y = item["value"]
            assert real_value in values
            assert y == 0
            assert item["name"] in dataset.all_samples


def test_get_option_violin_axis_has_no_static_data():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert all("data" not in x_axis for x_axis in option["xAxis"])
    assert all("data" not in y_axis for y_axis in option["yAxis"])


def test_violin_mark_count_is_total_scatter_points():
    plot = _make_violin_plot()
    dataset = plot.datasets[0]
    expected = sum(len(dataset.violin_value_by_sample_by_metric[m]) for m in dataset.metrics)
    assert echarts.violin.mark_count(dataset) == expected == 30  # 3 metrics * 10 samples


def test_serialize_violin_includes_violins_payload():
    plot = _make_violin_plot()
    result = echarts.serialize(plot)

    violins = result["datasets"][0]["violins"]
    dataset = plot.datasets[0]
    assert set(violins.keys()) == set(dataset.metrics)
    for metric in dataset.metrics:
        payload = violins[metric]
        assert len(payload["poly"]) == 2 * echarts.violin.N_BINS
        for x, y in payload["poly"]:
            assert isinstance(x, float)
            assert isinstance(y, float)
        header = dataset.header_by_metric[metric]
        values = echarts.violin._numeric_values(dataset, metric)
        assert payload["range"] == list(echarts.violin._metric_range(header, values))


def test_get_option_violin_inner_box_and_median():
    """POLISH.md #10b: each violin's `renderItem` also draws a Q1-Q3 box and median
    line, at the same normalized position `_row_quartiles` computes independently."""
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    dataset = plot.datasets[0]

    violin_series = [s for s in option["series"] if s["type"] == "custom" and "polygon" in _renderitem_body(s)]
    assert len(violin_series) == len(dataset.metrics)

    for metric, series_option in zip(dataset.metrics, violin_series):
        body = _renderitem_body(series_option)
        assert "'rect'" in body
        assert "'line'" in body

        values = echarts.violin._numeric_values(dataset, metric)
        q1, median, q3 = echarts.violin._row_quartiles(values)
        # The box/median coordinates are baked into the renderItem body as REAL x
        # literals (see `_violin_render_series`, no more per-row 0..1 normalization,
        # POLISH.md #6); check they actually appear in the generated source rather than
        # re-parsing the JS.
        for x in (q1, median, q3):
            assert json.dumps(x) in body


def test_row_quartiles_matches_box_quantile_method():
    """`_row_quartiles` reuses `box.py`'s quantile method directly (not a re-derivation),
    so this just pins the linear-interpolation contract for a simple known input."""
    q1, median, q3 = echarts.violin._row_quartiles([1.0, 2.0, 3.0, 4.0, 5.0])
    assert (q1, median, q3) == (2.0, 3.0, 4.0)


def test_metric_range_prefers_header_xaxis_range():
    """POLISH.md #10c: a column with a configured axis range (e.g. a 0-100 percentage)
    normalizes against that range, not the padded observed-value range, so the violin
    occupies the same fraction of its row that Plotly's does."""
    header = violin.ViolinColumn(
        title="Pct",
        description="",
        suffix="%",
        dmin=0,
        dmax=100,
        hidden=False,
        xaxis=violin.XAxis(range=[0.0, 100.0]),
        show_only_outliers=False,
        show_points=True,
    )
    assert echarts.violin._metric_range(header, [1.0, 2.0, 3.0]) == (0.0, 100.0)


def test_metric_range_falls_back_and_centers_degenerate_values():
    """No configured range: falls back to the padded observed range; when every value is
    identical, centers a small synthetic window on it instead of the old behavior of
    pinning it to the row's left edge (POLISH.md #10c/#10d)."""
    header = violin.ViolinColumn(
        title="Metric",
        description="",
        suffix="",
        dmin=None,
        dmax=None,
        hidden=False,
        xaxis=violin.XAxis(),
        show_only_outliers=False,
        show_points=True,
    )
    lo, hi = echarts.violin._metric_range(header, [5.0, 5.0, 5.0])
    assert lo < 5.0 < hi


def test_violin_polygon_degenerate_metric_is_flat():
    """POLISH.md #10d: a zero-variance metric must not blow kde()'s bandwidth floor
    into a needle-thin, full-height spike; the polygon should collapse to a flat
    (zero-height) line instead."""
    poly = echarts.violin._violin_polygon([5.0, 5.0, 5.0], lo=4.5, hi=5.5)
    ys = [y for _, y in poly]
    assert ys == [0.0, 0.0]


def test_violin_polygon_uses_real_x_not_normalized():
    """POLISH.md #6: each row draws on its own real-valued x-axis now, so the polygon's
    x coordinates are the metric's actual values, not fractions normalized to 0..1."""
    values = [10.0, 20.0, 30.0]
    lo, hi = echarts.violin._metric_range(
        violin.ViolinColumn(
            title="M",
            description="",
            suffix="",
            dmin=None,
            dmax=None,
            hidden=False,
            xaxis=violin.XAxis(),
            show_only_outliers=False,
            show_points=True,
        ),
        values,
    )
    poly = echarts.violin._violin_polygon(values, lo, hi)
    xs = [x for x, _ in poly]
    assert min(xs) >= lo
    assert max(xs) <= hi
    # Definitely not squeezed into 0..1 (the old normalized convention).
    assert max(xs) > 1.0


def test_row_geometry_rows_tile_without_overlap():
    """`_row_geometry` (PLOTLY-STYLE PER-ROW SUBPLOTS) must produce non-overlapping,
    monotonically increasing row slots that fit within the ideal 0-100% container."""
    n = 5
    prev_bottom = 0.0
    for row_idx in range(n):
        top_s, height_s = echarts.violin._row_geometry(row_idx, n)
        assert top_s.endswith("%")
        assert height_s.endswith("%")
        top, height = float(top_s[:-1]), float(height_s[:-1])
        assert top >= prev_bottom
        assert height > 0
        assert top + height <= 100.0
        prev_bottom = top  # next row's top must be >= this row's top (rows are ordered)


def test_grid_left_scales_with_title_length_and_is_clamped():
    """`_grid_left` (shared left inset so every row's x-axis starts at the same pixel,
    see the SSR-vs-interactive containLabel note on the function) grows with the
    longest title but stays within its documented floor/ceiling."""
    short = echarts.violin._grid_left(["A"])
    long = echarts.violin._grid_left(["A very much longer metric title than the others"])
    assert echarts.violin._MIN_GRID_LEFT <= short < long <= echarts.violin._MAX_GRID_LEFT
    assert echarts.violin._grid_left([]) == echarts.violin._MIN_GRID_LEFT


def test_get_option_violin_series_bound_to_own_row_axes():
    """Every series (violin/scatter) must carry `xAxisIndex`/`yAxisIndex` matching its
    row, so it draws in that row's own grid (PLOTLY-STYLE PER-ROW SUBPLOTS) rather than
    defaulting to row 0."""
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    dataset = plot.datasets[0]
    n = len(dataset.metrics)

    violin_series = [s for s in option["series"] if _is_kde_series(s)]
    scatter_series = [s for s in option["series"] if s["type"] == "scatter"]

    for row_idx, s in enumerate(violin_series):
        assert s["xAxisIndex"] == row_idx
        assert s["yAxisIndex"] == row_idx
    for row_idx, s in enumerate(scatter_series):
        assert s["xAxisIndex"] == row_idx
        assert s["yAxisIndex"] == row_idx
    assert len(violin_series) == len(scatter_series) == n


def test_get_option_violin_title_array_matches_rows():
    """ROW ALIGNMENT fix: the row title is a native `title` component entry (not a
    series), appended after the chart's own main title, one per visible metric, each
    fixed at the shared `grid_left` inset and vertically centered on its row."""
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    dataset = plot.datasets[0]
    n = len(dataset.metrics)

    titles = option["title"]
    assert isinstance(titles, list)
    assert len(titles) == n + 1  # main chart title, then one per row
    main_title, row_titles = titles[0], titles[1:]
    assert main_title.get("text")  # the chart's own "Test: Violin" title, untouched

    grid_lefts = {g["left"] for g in option["grid"]}
    assert len(grid_lefts) == 1  # every row's grid shares the same left (see the grid test)
    (grid_left,) = grid_lefts

    for row_idx, (metric, title_opt) in enumerate(zip(dataset.metrics, row_titles)):
        header = dataset.header_by_metric[metric]
        assert title_opt["text"] == echarts.violin._metric_title(header)
        assert title_opt["textAlign"] == "right"
        assert title_opt["textVerticalAlign"] == "middle"
        assert title_opt["padding"] == 0
        # Fixed pixel/percent position, not tied to any axis or series: same left for
        # every row (right-aligned just before that row's shared grid_left), vertically
        # centered on its own row.
        assert title_opt["left"] == grid_left - echarts.violin._TITLE_GUTTER_PAD
        assert title_opt["top"] == f"{echarts.violin._row_center_pct(row_idx, n):.4f}%"


############################################
# Task 3.2 (multiqc-echarts-exploration/BUILD_PLAN.md "test suite integration"):
# engine-parametrized tests, a box sort-by-median unit test, and end-to-end pipeline
# tests.


def _verify_rendered_for_engine(plot, plotting_engine: str) -> Dict[str, Any]:
    """
    Cloned subset of `tests/test_plots.py::_verify_rendered`, parametrized over both
    backends via the `plotting_engine` fixture (`tests/conftest.py`): render the plot and
    check the payload contract for whichever engine is active.
    """
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))
    dumped = report.plot_data[plot.anchor]
    if plotting_engine == "echarts":
        assert "echarts" in dumped
        payload = dumped["echarts"]
        assert "unsupported" not in payload
        assert payload["renderer"] in ("svg", "canvas")
        assert len(payload["datasets"]) == len(plot.datasets)
    else:
        assert "echarts" not in dumped
    return dumped


def test_engine_parametrized_bar(plotting_engine):
    _verify_rendered_for_engine(_make_bar_plot(), plotting_engine)


def test_engine_parametrized_line(plotting_engine):
    _verify_rendered_for_engine(_make_line_plot(), plotting_engine)


def test_engine_parametrized_scatter(plotting_engine):
    _verify_rendered_for_engine(_make_scatter_plot(), plotting_engine)


def test_engine_parametrized_heatmap(plotting_engine):
    _verify_rendered_for_engine(_make_heatmap_plot(), plotting_engine)


def test_engine_parametrized_box(plotting_engine):
    _verify_rendered_for_engine(_make_box_plot(), plotting_engine)


def test_engine_parametrized_violin(plotting_engine):
    _verify_rendered_for_engine(_make_violin_plot(), plotting_engine)


def test_box_sort_by_median_produces_median_sorted_series():
    """
    D3 (`multiqc-echarts-exploration/PARITY.md`): the sort-by-median toggle is wired up
    only on the JS side (the inherited, engine-neutral `prepData()` swaps in
    `dataset.data_sorted`/`samples_sorted` when the sort switch is active,
    `templates/default/src/js/plots-echarts/box.js`); no shipped module sets `sort_by_median=True`
    and nothing exercises that sorted path. `Dataset.create` (`multiqc/plots/box.py`)
    always keeps `dataset.data`/`samples` as the alphabetical order and stashes the
    ascending-by-median order separately in `data_sorted`/`samples_sorted`. Simulate
    exactly what the JS toggle does (swap those in as the active dataset fields) and
    drive the result through the real `echarts.get_option()` used by the SSR/flat export
    path, proving the ECharts box builder (which only ever reads `dataset.data`/
    `samples`) is order-agnostic and renders a genuinely median-sorted series.
    """
    data = {
        "SampleB": [1.0, 2.0, 3.0],  # median 2
        "SampleA": [10.0, 20.0, 30.0],  # median 20
        "SampleC": [100.0, 200.0, 300.0],  # median 200
    }
    plot = box.plot(
        data,
        box.BoxPlotConfig(id="boxplot_sorted", title="Test: Box Sorted", sort_by_median=True),
    )
    assert isinstance(plot, box.BoxPlot)
    dataset = plot.datasets[0]
    assert dataset.samples_sorted == ["SampleB", "SampleA", "SampleC"]
    assert dataset.data_sorted is not None

    sorted_dataset = box.Dataset(**{**dataset.__dict__, "samples": dataset.samples_sorted, "data": dataset.data_sorted})
    plot.datasets = [sorted_dataset]

    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["yAxis"]["data"] == ["SampleB", "SampleA", "SampleC"]
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    medians = [entry["value"][2] for entry in boxplot_series["data"]]
    assert medians == [2.0, 20.0, 200.0]
    assert medians == sorted(medians)


############################################
# End-to-end: the full `multiqc.run()` pipeline against a single, self-contained module
# data file (fastp's SAMPLE.json; no directory scan needed) so this stays fast. fastp
# emits bar, line, and violin plots, so this also touches most of the ported plot types.


def test_echarts_template_end_to_end_interactive(data_dir, tmp_path):
    result = multiqc.run(
        data_dir / "modules" / "fastp" / "SAMPLE.json",
        cfg=ClConfig(run_modules=["fastp"], template="default", output_dir=tmp_path),
        return_html=True,
    )

    assert result.sys_exit_code == 0
    assert result.html_content is not None
    # The default template's footer credits Apache ECharts, proving the ECharts bundle
    # was actually used (default renders ECharts).
    assert "Apache ECharts" in result.html_content

    assert len(report.plot_data) > 0
    for dumped in report.plot_data.values():
        assert "echarts" in dumped
        assert "unsupported" not in dumped["echarts"]


def test_echarts_template_end_to_end_flat(data_dir, tmp_path):
    result = multiqc.run(
        data_dir / "modules" / "fastp" / "SAMPLE.json",
        cfg=ClConfig(run_modules=["fastp"], template="default", plots_force_flat=True, output_dir=tmp_path),
        return_html=True,
    )

    assert result.sys_exit_code == 0
    assert result.html_content is not None
    # Static export path (SSR + resvg PNG rasterization): flat plots are embedded as
    # base64 PNGs, mirroring the Plotly flat-plot contract.
    assert "data:image/png;base64," in result.html_content
