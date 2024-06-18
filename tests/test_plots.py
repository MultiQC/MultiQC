from multiqc import report
from multiqc.plots import bargraph, linegraph, table, violin, heatmap, scatter, box


def test_barplot():
    plot_id = "test_barplot"

    plot = bargraph.plot(
        {
            "Sample0": {},
            "Sample1": {"Cat1": 1},
            "Sample2": {"Cat1": 1, "Cat2": 1},
            "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
        },
        ["Cat1", "Cat2"],
        {
            "id": plot_id,
            "title": "Test: Bar Graph",
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_plot_no_matching_data():
    plot_id = "test_plot_no_matching_data"

    plot = bargraph.plot(
        {
            "Sample1": {"Cat0": 1, "Cat1": 1},
        },
        ["Cat2", "Cat3"],
        {
            "id": plot_id,
            "title": "Test: Bar Graph",
        },
    )

    assert isinstance(plot, str)


def test_linegraph():
    plot_id = "test_linegraph"

    plot = linegraph.plot(
        {
            "Sample0": {},
            "Sample1": {0: 1, 1: 1},
            "Sample2": {0: 1, 1: 1, 2: 1},
        },
        {
            "id": plot_id,
            "title": "Test: Line Graph",
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_multiple_datasets():
    plot_id = "test_multiple_datasets"

    plot = linegraph.plot(
        [
            {
                "Sample1": {0: 1, 1: 1},
            },
            {
                "Sample1": {0: 2, 1: 2},
            },
        ],
        {
            "id": plot_id,
            "title": "Test: Line Graph",
            "data_labels": ["Dataset1", "Dataset2"],
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


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
        pconfig={
            "id": plot_id,
            "title": "Test: Table",
        },
    )

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
        pconfig={
            "id": plot_id,
            "title": "Test: Violin Plot",
        },
    )

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
        pconfig={
            "id": plot_id,
            "title": "Test: Heatmap",
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_scatter():
    plot_id = "test_scatter"

    plot = scatter.plot(
        {
            "Sample1": [{"x": 1, "y": 2}],
        },
        {
            "id": plot_id,
            "title": "Test: Scatter Plot",
            "xlab": "X",
            "ylab": "Y",
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def text_box():
    plot_id = "test_box"

    plot = box.plot(
        {
            "Sample1": [1, 2, 3],
            "Sample2": [4, 5, 6],
        },
        {
            "id": plot_id,
            "title": "Test: Box Plot",
        },
    )

    report.reset()
    plot.add_to_report()
    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data


def test_flat_plot():
    plot_id = "test_flat_plot"

    plot = linegraph.plot(
        {
            "Sample1": {0: 1, 1: 1},
        },
        {
            "id": plot_id,
            "title": "Test: Flat Line Graph",
        },
    )

    plot.flat = True
    report.reset()
    html = plot.add_to_report()
    assert len(report.plot_data) == 0
    assert html is not None
    assert f"""<div class="mqc_mplplot" id="{plot_id}"><img src="data:image/png;base64""" in html
