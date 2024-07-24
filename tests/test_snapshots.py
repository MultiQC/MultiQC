import pytest

from multiqc import report, BaseMultiqcModule, config
from multiqc.core.write_results import write_results
from multiqc.plots import table


@pytest.fixture
def custom_module():
    module = BaseMultiqcModule(name="my-module", anchor="custom_data")
    module.add_section(
        name="Custom Section",
        description="Custom description",
        helptext="Custom help",
        plot=table.plot(
            data={"sample1": {"x": 1, "y": 2}, "sample2": {"x": 3, "y": 4}},
            headers={"x": {"title": "Metric X"}},
            pconfig=table.TableConfig(id="custom_table", title="Custom table"),
        ),
    )
    module.general_stats_addcols(
        data={"sample1": {"x": 42}},
        headers={"x": {"title": "Custom title"}},
    )
    report.modules = [module]


@pytest.fixture
def browser():
    from selenium import webdriver

    driver = webdriver.Chrome()  # You can use other drivers like Firefox, Edge, etc.
    yield driver
    driver.quit()


def test_custom_module(custom_module, tmp_path, snapshot, browser):
    """
    Verify snapshot option: HTML is written to disk, JSON is not
    """
    config.version = "0.0.0"
    config.creation_date = "2021-01-01"
    config.make_data_dir = False
    config.output_fn = tmp_path / "report.html"
    config.development = True

    # Renders plots and writes HTML, but does not re-init config like multiqc.write_report()
    write_results()

    screenshot_path = tmp_path / "screenshot.png"
    print(screenshot_path)

    # Load the report in the browser and take a screenshot
    browser.get(f"file://{config.output_fn}")
    browser.save_screenshot(screenshot_path)
    assert screenshot_path.is_file()

    snapshot.assert_match(screenshot_path.read_bytes(), "screenshot.png")


# TODO: test custom content and other modules and plot types
