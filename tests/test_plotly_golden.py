"""
Golden regression guard for the Plotly rendering path.

Pins, per plot type, (a) the interactive browser payload `plot.model_dump()["layout"]`
and (b) the full static `get_figure(0).to_plotly_json()` (layout + traces). This is the
safety net for migrating `Plot.layout` off `go.Layout` onto the neutral IR: the layout the
browser Plotly renderer and the static export receive must stay byte-identical through the
refactor. Regenerate deliberately with `UPDATE_GOLDEN=1 pytest tests/test_plotly_golden.py`
(and eyeball the diff) whenever a change is intended.
"""

import json
import os
from pathlib import Path

import pytest

from multiqc.plots import bargraph, box, heatmap, linegraph, scatter, seqcontent, table, violin

GOLDEN_DIR = Path(__file__).parent / "data" / "plotly_golden"


def _fastqc_shaped_data():
    return {
        "sample_2": {
            "1": {"a": 25.0, "c": 25.0, "g": 25.0, "t": 25.0},
            "2-3": {"a": 20.0, "c": 30.0, "g": 30.0, "t": 20.0},
        },
        "sample_1": {"1": {"a": 10.0, "c": 20.0, "g": 30.0, "t": 40.0}},
    }


def _build(name):
    if name == "bar":
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
    if name == "line":
        return linegraph.plot(
            {"Sample0": {}, "Sample1": {0: 1, 1: 1}, "Sample2": {0: 1, 1: 1, 2: 1}},
            linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
        )
    if name == "violin":
        return violin.plot(
            data={"sample1": {"x": 1, "y": 2}, "sample2": {"x": 3, "y": 4}},
            headers={"x": {"title": "Metric X"}},
            pconfig=table.TableConfig(id="violin", title="Violin"),
        )
    if name == "heatmap":
        return heatmap.plot(
            data=[[1, 2], [3, 4]],
            xcats=["Cat1", "Cat2"],
            ycats=["Sample1", "Sample2"],
            pconfig=heatmap.HeatmapConfig(id="heatmap", title="Heatmap"),
        )
    if name == "scatter":
        return scatter.plot(
            {"Sample1": [{"x": 1, "y": 2}]},
            scatter.ScatterConfig(id="scatter", title="Scatter", xlab="X", ylab="Y"),
        )
    if name == "box":
        return box.plot(
            {"Sample1": [1, 2, 3], "Sample2": [4, 5, 6]},
            box.BoxPlotConfig(id="box", title="Box"),
        )
    if name == "seqcontent":
        return seqcontent.plot(_fastqc_shaped_data(), {"id": "seqcontent_test", "title": "Test seqcontent"})
    raise ValueError(name)


NAMES = ["bar", "line", "violin", "heatmap", "scatter", "box", "seqcontent"]


def _normalize(obj):
    """Round-trip through JSON (default=str) so numpy/tuples/ordering are comparable."""
    return json.loads(json.dumps(obj, default=str, sort_keys=True))


@pytest.mark.parametrize("name", NAMES)
def test_plotly_layout_golden(name):
    plot = _build(name)
    got = _normalize(
        {
            "interactive_layout": plot.model_dump(warnings=False)["layout"],
            "figure": plot.get_figure(0).to_plotly_json(),
        }
    )
    golden_path = GOLDEN_DIR / f"{name}.json"
    if os.environ.get("UPDATE_GOLDEN"):
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(json.dumps(got, indent=2, sort_keys=True))
        pytest.skip(f"wrote golden {golden_path.name}")
    golden = json.loads(golden_path.read_text())
    assert got == golden
