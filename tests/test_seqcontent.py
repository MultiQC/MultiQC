import plotly.graph_objects as go
import pytest

from multiqc.core.special_case_modules.load_multiqc_data import create_plot_input_data_only
from multiqc.plots import seqcontent
from multiqc.plots.seqcontent import (
    SeqContentBin,
    SeqContentNormalizedInputData,
    SeqContentPlot,
    _parse_bin_label,
    bin_rgb,
)
from multiqc.types import PlotType


def _fastqc_shaped_data():
    """Two samples, shaped like fastqc's data_by_sample: sample -> bin label -> a/c/g/t."""
    return {
        "sample_2": {
            "1": {"a": 25.0, "c": 25.0, "g": 25.0, "t": 25.0},
            "2-3": {"a": 20.0, "c": 30.0, "g": 30.0, "t": 20.0},
        },
        "sample_1": {
            "1": {"a": 10.0, "c": 20.0, "g": 30.0, "t": 40.0},
        },
    }


def test_basic_creation_from_fastqc_shaped_dict():
    plot = seqcontent.plot(
        _fastqc_shaped_data(),
        {"id": "seqcontent_test", "title": "Test seqcontent"},
    )
    assert isinstance(plot, SeqContentPlot)
    assert len(plot.datasets) == 1
    ds = plot.datasets[0]
    # natsorted sample order
    assert ds.samples == ["sample_1", "sample_2"]
    assert len(ds.rows) == 2
    assert ds.rows[1][1].label == "2-3"
    assert ds.rows[1][1].start == 2
    assert ds.rows[1][1].end == 3
    assert ds.max_bp == 3


def test_rgb_golden_t_100():
    b = SeqContentBin(label="1", start=1, end=1, a=0.0, c=0.0, g=0.0, t=100.0)
    r, g, bl = bin_rgb(b)
    assert (r, g, bl) == (255, 0, 0)


def test_rgb_golden_even_split():
    b = SeqContentBin(label="1", start=1, end=1, a=25.0, c=25.0, g=25.0, t=25.0)
    r, g, bl = bin_rgb(b)
    assert (r, g, bl) == (64, 64, 64)


def test_old_count_based_renorm():
    """Very old FastQC gives raw counts instead of percentages: t=300,a=100,c=100,g=100
    (t > 100) triggers renormalization to percentages of the sum (600):
    t=50%, a=c=g=16.67%."""
    data_by_sample = {"sample1": {"1": {"a": 100.0, "c": 100.0, "g": 100.0, "t": 300.0}}}
    inputs = SeqContentNormalizedInputData.create(data_by_sample, {"id": "x", "title": "x"})
    b = inputs.data["sample1"][0]
    assert b.t == 50.0
    assert b.a == pytest.approx(16.67)
    assert b.c == pytest.approx(16.67)
    assert b.g == pytest.approx(16.67)
    r, g, bl = bin_rgb(b)
    assert r == 128
    assert g == round(16.67 / 100 * 255)
    assert bl == round(16.67 / 100 * 255)


def test_nan_and_missing_base_become_zero():
    data_by_sample = {
        "sample1": {
            "1": {"a": float("nan"), "c": 10.0, "g": 10.0, "t": 10.0},
            "2": {"c": 10.0, "g": 10.0, "t": 10.0},  # "a" missing entirely
        }
    }
    inputs = SeqContentNormalizedInputData.create(data_by_sample, {"id": "x", "title": "x"})
    bins = inputs.data["sample1"]
    assert bins[0].a == 0.0
    assert bins[1].a == 0.0


def test_bin_label_parsing():
    assert _parse_bin_label("7") == (7, 7)
    assert _parse_bin_label("2-3") == (2, 3)
    assert _parse_bin_label("10-14") == (10, 14)


def test_bin_label_parsing_crashes_on_garbage():
    with pytest.raises(ValueError):
        _parse_bin_label("not-a-bin")


def test_to_df_from_df_round_trip():
    inputs = SeqContentNormalizedInputData.create(_fastqc_shaped_data(), {"id": "roundtrip", "title": "Roundtrip"})
    df = inputs.to_df()
    assert not df.is_empty()

    reloaded = SeqContentNormalizedInputData.from_df(df, inputs.pconfig, inputs.anchor)
    assert reloaded.data.keys() == inputs.data.keys()
    for s_name in inputs.data:
        orig_bins = inputs.data[s_name]
        new_bins = reloaded.data[s_name]
        assert [b.label for b in orig_bins] == [b.label for b in new_bins]
        assert [b.t for b in orig_bins] == [b.t for b in new_bins]


def test_merge_sample_union():
    old = SeqContentNormalizedInputData.create(
        {"sampleA": {"1": {"a": 10, "c": 10, "g": 10, "t": 70}}},
        {"id": "merge_test", "title": "Merge"},
    )
    new = SeqContentNormalizedInputData.create(
        {"sampleB": {"1": {"a": 20, "c": 20, "g": 20, "t": 40}}},
        {"id": "merge_test", "title": "Merge"},
    )
    merged = SeqContentNormalizedInputData.merge(old, new)
    assert set(merged.data.keys()) == {"sampleA", "sampleB"}

    # New run wins per-sample: overlapping sample uses the new data.
    updated = SeqContentNormalizedInputData.create(
        {"sampleA": {"1": {"a": 1, "c": 1, "g": 1, "t": 97}}},
        {"id": "merge_test", "title": "Merge"},
    )
    merged2 = SeqContentNormalizedInputData.merge(merged, updated)
    assert merged2.data["sampleA"][0].t == 97.0
    assert "sampleB" in merged2.data


def test_model_dump_json_round_trip_through_create_plot_input_data_only():
    inputs = SeqContentNormalizedInputData.create(
        _fastqc_shaped_data(), {"id": "json_roundtrip", "title": "JSON roundtrip"}
    )
    # Mirrors save_to_parquet(): plot_type serializes to its enum value in JSON mode.
    dumped = inputs.model_dump(mode="json", exclude_none=True)

    reloaded = create_plot_input_data_only(dumped)
    assert isinstance(reloaded, SeqContentNormalizedInputData)
    assert reloaded.plot_type == PlotType.SEQCONTENT
    assert reloaded.data.keys() == inputs.data.keys()


def test_create_figure_returns_single_image_trace():
    plot = seqcontent.plot(
        _fastqc_shaped_data(),
        {"id": "figure_test", "title": "Figure test"},
    )
    assert isinstance(plot, SeqContentPlot)
    fig = plot.get_figure(0, flat=True)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1
    trace = fig.data[0]
    assert isinstance(trace, go.Image)
    n_samples = len(plot.datasets[0].samples)
    max_bp = plot.datasets[0].max_bp
    assert trace.z.shape == (n_samples, max_bp, 3)


def test_static_export_sets_yaxis_scaleratio_to_stretch_image():
    """
    go.Image traces force an equal-aspect y axis (yaxis.scaleanchor) that kaleido
    cannot clear, collapsing the heatmap to a thin strip unless yaxis.scaleratio is
    set (see seqcontent.py::_static_scaleratio, twin of
    seqcontent.js::fixAspectRatio). The flat/static figure must set it to a
    positive value for a multi-sample dataset.
    """
    plot = seqcontent.plot(
        _fastqc_shaped_data(),
        {"id": "scaleratio_test", "title": "Scaleratio test"},
    )
    assert isinstance(plot, SeqContentPlot)
    fig = plot.get_figure(0, flat=True)
    assert fig.layout.yaxis.scaleratio is not None
    assert fig.layout.yaxis.scaleratio > 0
