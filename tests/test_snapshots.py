import asyncio
import time
from typing import Dict, Union

import pytest

from multiqc import report, BaseMultiqcModule, config, Plot
from multiqc.core.write_results import write_results
from multiqc.plots import bargraph, linegraph, box, table, violin, heatmap, scatter

from pyppeteer import launch  # type: ignore


def _plots() -> Dict[str, Plot]:
    """
    Custom plots that we want to screenshot-test. One for each plot type.
    """

    return {
        str(p.plot_type): p
        for p in [
            bargraph.plot(
                {
                    "Sample0": {},
                    "Sample1": {"Cat1": 1},
                    "Sample2": {"Cat1": 1, "Cat2": 1},
                    "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
                },
                ["Cat1", "Cat2"],
                bargraph.BarPlotConfig(id="bargraph", title="Test: Bar Graph"),
            ),
            linegraph.plot(
                {
                    "Sample0": {},
                    "Sample1": {0: 1, 1: 1},
                    "Sample2": {0: 1, 1: 1, 2: 1},
                },
                linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
            ),
            table.plot(
                data={
                    "sample1": {"x": 1, "y": 2},
                    "sample2": {"x": 3, "y": 4},
                },
                headers={"x": {"title": "Metric X"}},
                pconfig=table.TableConfig(id="table", title="Table"),
            ),
            violin.plot(
                data={"sample1": {"x": 1, "y": 2}, "sample2": {"x": 3, "y": 4}},
                headers={"x": {"title": "Metric X"}},
                pconfig=table.TableConfig(id="violin", title="Violin"),
            ),
            heatmap.plot(
                data=[[1, 2], [3, 4]],
                xcats=["Cat1", "Cat2"],
                ycats=["Sample1", "Sample2"],
                pconfig=heatmap.HeatmapConfig(id="heatmap", title="Heatmap"),
            ),
            scatter.plot(
                {"Sample1": [{"x": 1, "y": 2}]},
                scatter.ScatterConfig(id="scatter", title="Scatter", xlab="X", ylab="Y"),
            ),
            box.plot(
                {"Sample1": [1, 2, 3], "Sample2": [4, 5, 6]},
                box.BoxPlotConfig(id="box", title="Box"),
            ),
            None,  # general stats
        ]
        if isinstance(p, Plot)
    }


def _add_module(plot: Union[Plot, None]):
    """
    Add module to report. Should be called inside the test function as it modifies the `report` object
    (in two places: the `module.add_section` call and in `report.modules = [module]`)
    """
    module = BaseMultiqcModule(name="Custom Module", anchor="custom_module")
    if plot is None:
        module.general_stats_addcols(
            data={"sample1": {"x": 42}},
            headers={"x": {"title": "Custom title"}},
        )
    else:
        module.add_section(
            name=f"Custom Section {plot.id}",
            description=f"Custom description for {plot.id}",
            helptext=f"Custom help for {plot.id}",
            plot=plot,
        )
    report.modules = [module]


async def get_full_page_screenshot(html_path, out_png_path):
    browser = await launch(headless=True, args=["--no-sandbox"])
    page = await browser.newPage()
    await page.goto(f"file://{html_path}", waitUntil="networkidle0")
    await page.setViewport({"width": 1080, "height": 1080})
    await page.screenshot({"path": out_png_path, "fullPage": True})
    await browser.close()


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "name,plot",
    list(_plots().items()) + [("general_stats", None)],
)
async def test_custom_module_with_plot(tmp_path, snapshot, name, plot: Union[Plot, None]):
    """
    Verify snapshot option: HTML is written to disk, JSON is not
    """
    config.version = "0.0.0"
    config.creation_date = "2021-01-01"
    config.make_data_dir = False
    config.output_dir = tmp_path
    config.development = True  # no need to compress and embed JS

    _add_module(plot)

    # Renders plots and writes HTML, but does not re-init config like multiqc.write_report()
    write_results()

    screenshot_path = tmp_path / "screenshot.png"
    await get_full_page_screenshot(tmp_path / "multiqc_report.html", screenshot_path)
    assert screenshot_path.is_file()

    print(screenshot_path)
    assert screenshot_path.read_bytes() == snapshot, screenshot_path
