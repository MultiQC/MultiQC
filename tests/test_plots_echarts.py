"""
Serialization contract tests for the ECharts backend (`multiqc/plots/echarts/`).

Uses the same bar plot as `tests/test_plots.py::test_barplot`. The autouse `reset`
fixture in `tests/conftest.py` resets `config`/`report` after every test.
"""

import math

from multiqc import config, report
from multiqc.plots import bargraph, echarts, linegraph
from multiqc.types import Anchor


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


def _make_line_plot():
    return linegraph.plot(
        {
            "Sample1": {0: 1, 1: 1},
            "Sample2": {0: 1, 1: 1, 2: 1},
        },
        linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
    )


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
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["title"]["text"] == "Foo Bar"
    assert skeleton["title"]["subtext"] == "42 things"


def test_interactive_plot_unsupported_plot_type_does_not_raise():
    config.plotting_engine = "echarts"
    plot = _make_line_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert dumped["echarts"] == {"unsupported": "x/y line"}


def test_interactive_plot_omits_echarts_key_by_default():
    assert config.plotting_engine == "plotly"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" not in dumped


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


def test_get_option_is_log_sets_switch_controlled_axis_type():
    plot = _make_bar_plot()
    assert plot.axis_controlled_by_switches == ["xaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["xAxis"]["type"] == "log"
