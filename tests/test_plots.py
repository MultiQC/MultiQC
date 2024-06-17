from multiqc import report
from multiqc.plots import bargraph


def test_barplot():
    """
    Verify that all modules do at least something
    """
    report.reset()

    plot_id = "test_barplot"

    plot = bargraph.plot(
        {
            "Sample1": {"Cat1": 10, "Cat2": 20},
            "Sample2": {"Cat1": 5, "Cat2": 15, "Cat3": 25},
        },
        ["Cat1", "Cat2", "Cat3"],
        {
            "id": plot_id,
            "title": "Test: Bar Graph",
        },
    )

    plot.add_to_report()

    assert len(report.plot_data) == 1
    assert plot_id in report.plot_data
