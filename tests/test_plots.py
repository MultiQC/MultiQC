import pytest

from multiqc import report, Plot
from multiqc.core.exceptions import RunError
from multiqc.plots import bargraph, linegraph, table, violin, heatmap, scatter, box


def test_barplot():
    """
    Normal barplot call. Single datasets, a list of cats, some cats are excluded.
    """
    plot_id = "test_barplot"

    plot = bargraph.plot(
        {
            "Sample0": {},
            "Sample1": {"Cat1": 1},
            "Sample2": {"Cat1": 1, "Cat2": 1},
            "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
        },
        ["Cat1", "Cat2"],
        {"id": plot_id, "title": "Test: Bar Graph"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_bar_plot_no_matching_cats():
    """
    None of the cats are matching thd data, so shouldn't produce a plot
    """
    plot_id = "test_bar_plot_no_matching_cats"

    plot = bargraph.plot(
        {"Sample1": {"Cat0": 1, "Cat1": 1}},
        ["Cat2", "Cat3"],
        {"id": plot_id, "title": "Test: Bar Graph"},
    )
    # Will return a warning message html instead of a plot:
    assert isinstance(plot, str)


def test_bar_plot_cats_dicts():
    """
    Advanced cats spec - dict with cat properties instead of a simple list
    """
    plot_id = "test_bar_plot_cats_dicts"

    plot = bargraph.plot(
        {"Sample1": {"Cat1": 1}},
        {"Cat1": {"name": "My category"}},
        {"id": plot_id, "title": "Test: Bar Graph"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert report.plot_data[plot_id]["datasets"][0]["cats"][0]["name"] == "My category"


def test_bar_plot_cats_dicts_with_typo():
    """
    A typo in the cat properties dict - will fill in the name from the dict key
    """
    plot_id = "test_bar_plot_cats_dicts_with_typo"

    plot = bargraph.plot(
        {"Sample1": {"Cat1": 2}},
        {"Cat1": {"name_with_typo": "My category"}},
        {"id": plot_id, "title": "Test: Bar Graph"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert report.plot_data[plot_id]["datasets"][0]["cats"][0]["name"] == "Cat1"


def test_bar_plot_cats_mismatch_cats_and_ds_count():
    """
    Multiple datasets, but the list lengths are not matching between data and cats. Should throw an error
    """
    plot_id = "test_bar_plot_cats_mismatch_cats_and_ds_count"

    with pytest.raises(RunError):
        bargraph.plot(
            [{"Sample1": {"Cat1": 2}}],
            [{"Cat1": {"name": "My category"}}, {"Cat1": {"name": "My category"}}],
            {"id": plot_id, "title": "Test: Bar Graph"},
        )


def test_bar_plot_fill_cats():
    """
    Multiple datasets, but only one dict of cats - in this case, should copy these cats between datasets
    """
    plot_id = "test_bar_plot_fill_cats"

    plot = bargraph.plot(
        [{"Sample1": {"Cat1": 2, "Cat2": 2}}, {"Sample1": {"Cat1": 2, "Cat3": 2}}],
        {"Cat1": {"name": "My category"}},
        {"id": plot_id, "title": "Test: Bar Graph"},
    )
    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data[plot_id]["datasets"]) == 2
    assert len(report.plot_data[plot_id]["datasets"][0]["cats"]) == 1
    assert len(report.plot_data[plot_id]["datasets"][1]["cats"]) == 1
    assert report.plot_data[plot_id]["datasets"][0]["cats"][0]["name"] == "My category"
    assert report.plot_data[plot_id]["datasets"][1]["cats"][0]["name"] == "My category"


def test_bar_plot_no_cats():
    """
    Cats parameter is missing, will fill in one from the data
    """

    plot_id = "test_bar_plot_no_cats"

    plot = bargraph.plot(
        {
            "Sample1": {"Cat1": 2, "Cat2": 2},
            "Sample2": {"Cat1": 1, "Cat3": 1},
        },
        pconfig={"id": plot_id, "title": "Test: Bar Graph"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data[plot_id]["datasets"][0]["cats"]) == 3


def test_linegraph():
    plot_id = "test_linegraph"

    dataset = {
        "Sample0": {},
        "Sample1": {0: 1, 1: 1},
        "Sample2": {0: 1, 1: 1, 2: 1},
    }

    plot = linegraph.plot(
        dataset,
        {"id": plot_id, "title": "Test: Line Graph"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data

    for in_series, out_series in zip(dataset.values(), report.plot_data[plot_id]["datasets"][0]["lines"]):
        assert len(in_series) == len(out_series["pairs"])


def test_linegraph_smooth():
    plot_id = "test_linegraph_smooth"
    SMOOTH_TO = 2
    dataset = {"Smoothed": {0: 1, 1: 1, 2: 1}, "Unsmoothed": {0: 1, 1: 1}}
    plot = linegraph.plot(
        dataset,
        {"id": plot_id, "title": "Test: Line Graph", "smooth_points": SMOOTH_TO},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    for in_series, out_series in zip(dataset.values(), report.plot_data[plot_id]["datasets"][0]["lines"]):
        assert min(len(in_series), SMOOTH_TO) == len(out_series["pairs"])


def test_linegraph_multiple_datasets():
    plot_id = "test_linegraph_multiple_datasets"

    plot = linegraph.plot(
        [{"Sample1": {0: 1, 1: 1}}, {"Sample1": {0: 2, 1: 2}}],
        {
            "id": plot_id,
            "title": "Test: Line Graph",
            "data_labels": ["Dataset1", "Dataset2"],
        },
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data[plot_id]["datasets"]) == 2


def test_table():
    plot_id = "test_table"

    plot = table.plot(
        {
            "Sample1": {"Metric1": 1, "Metric2": 2},
            "Sample2": {"Metric1": 3, "Metric2": 4},
        },
        headers={
            "Metric1": {"title": "Metric 1"},
            "Metric2": {"title": "Metric 2"},
        },
        pconfig={"id": plot_id, "title": "Test: Table"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_violin():
    plot_id = "test_violin"

    plot = violin.plot(
        {
            "Sample1": {"Metric1": 1, "Metric2": 2},
            "Sample2": {"Metric1": 3, "Metric2": 4},
        },
        headers={
            "Metric1": {"title": "Metric 1"},
            "Metric2": {"title": "Metric 2"},
        },
        pconfig={"id": plot_id, "title": "Test: Violin Plot"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_heatmap():
    plot_id = "test_heatmap"

    plot = heatmap.plot(
        data=[[1, 2], [3, 4]],
        xcats=["Cat1", "Cat2"],
        ycats=["Sample1", "Sample2"],
        pconfig={"id": plot_id, "title": "Test: Heatmap"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_scatter():
    plot_id = "test_scatter"

    plot = scatter.plot(
        {"Sample1": [{"x": 1, "y": 2}]},
        {
            "id": plot_id,
            "title": "Test: Scatter Plot",
            "xlab": "X",
            "ylab": "Y",
        },
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_box_plot():
    plot_id = "test_box_plot"

    plot = box.plot(
        {"Sample1": [1, 2, 3], "Sample2": [4, 5, 6]},
        {"id": plot_id, "title": "Test: Box Plot"},
    )

    assert isinstance(plot, Plot)
    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_flat_plot():
    plot_id = "test_flat_plot"

    plot = linegraph.plot(
        {"Sample1": {0: 1, 1: 1}},
        {"id": plot_id, "title": "Test: Flat Line Graph"},
    )

    assert isinstance(plot, Plot)
    plot.flat = True
    report.reset()
    html = plot.add_to_report()
    assert len(report.plot_data) == 0
    assert html is not None
    assert f"""<div class="mqc_mplplot" id="{plot_id}"><img src="data:image/png;base64""" in html
