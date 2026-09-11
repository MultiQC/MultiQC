import json
from pathlib import Path

import pytest

from multiqc import config, report, validation
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.nexons import MultiqcModule


@pytest.fixture(autouse=True)
def reset():
    report.reset()
    config.reset()
    validation.reset()
    yield
    report.reset()
    config.reset()
    validation.reset()


@pytest.fixture
def stats():
    return {
        "file": "/input/sample.bam",
        "outcomes": {
            "Total_Reads": 100,
            "No_Alignment": 10,
            "Primary_Alignment": 80,
            "Secondary_Alignment": 10,
            "No_Gene": 5,
            "No_Hit": 0,
            "Gene": 70,
            "Multi_Gene": 5,
            "Partial": 30,
            "Unique": 20,
            "Same_Strand_Hit": 60,
            "Opposing_Strand_Hit": 15,
        },
        "read_lengths": [[0, 10], [100, 70], [200, 20]],
        "coverage": [10] * 101,
        "inner_flex": {"2": 3, "-2": 1, "0": 96},
        "end_flex": {"-500": 10, "0": 80, "500": 10},
    }


def run_module(tmp_path: Path, stats):
    path = tmp_path / "output_nexons_stats.txt"
    path.write_text(json.dumps(stats))
    report.analysis_files = [path]
    report.search_files(["nexons"])
    return MultiqcModule()


def test_full_report(tmp_path, stats):
    module = run_module(tmp_path, stats)
    assert module.nexons_data["sample"] == stats
    assert len(module.sections) == 12
    assert module.percentages["sample"]["Gene"] == 70
    assert module.percentages["sample"]["Unique"] == 20
    # Directionality is divided by ALL reads, not just strand matches.
    assert module.percentages["sample"]["Same_Strand_Hit"] == 60
    assert module.percentages["sample"]["Opposing_Strand_Hit"] == 15
    assert not validation._warnings_by_cfg_path
    assert not validation._errors_by_cfg_path


def test_zero_counts(tmp_path, stats):
    stats["outcomes"] = dict.fromkeys(stats["outcomes"], 0)
    stats.update(read_lengths=[], coverage=[], inner_flex={}, end_flex={})
    module = run_module(tmp_path, stats)
    assert set(module.percentages["sample"].values()) == {0}
    assert len(module.sections) == 12
    assert module.sections[-1].alerts[0].affected_samples == ["sample"]


def test_ignore_samples(tmp_path, stats):
    config.sample_names_ignore = ["sample"]
    with pytest.raises(ModuleNoSamplesFound):
        run_module(tmp_path, stats)


def test_filename_sample_name(tmp_path, stats):
    config.use_filename_as_sample_name = True
    assert "output" in run_module(tmp_path, stats).nexons_data


def test_no_files():
    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()


@pytest.mark.parametrize("contents", ["", "not json", "[]", "{}"])
def test_invalid_json_or_structure(contents):
    with pytest.raises((ValueError, KeyError)):
        MultiqcModule.parse_stats(contents)


def test_missing_required_metric(stats):
    del stats["outcomes"]["Gene"]
    with pytest.raises(KeyError):
        MultiqcModule.parse_stats(json.dumps(stats))


@pytest.mark.parametrize("count", [-1, 1.5, "10", True, float("nan")])
def test_invalid_count(stats, count):
    stats["outcomes"]["Gene"] = count
    with pytest.raises(ValueError):
        MultiqcModule.parse_stats(json.dumps(stats))


def test_multiple_samples(tmp_path, stats):
    for i in range(2):
        stats["file"] = f"sample{i}.bam"
        (tmp_path / f"sample{i}_nexons_stats.txt").write_text(json.dumps(stats))
    report.analysis_files = [tmp_path]
    report.search_files(["nexons"])
    assert set(MultiqcModule().nexons_data) == {"sample0", "sample1"}


def test_plotted_values_preserve_bins_and_denominators(tmp_path, stats):
    stats["outcomes"]["Total_Reads"] = 200
    stats["end_flex"] = {str(x): x + 500 for x in range(-500, 501)}
    run_module(tmp_path, stats)
    strand = report.plot_by_id["nexons_directionality_plot"].datasets
    assert strand[0].cats[0].data == [30]
    assert strand[1].cats[0].data == [60]
    end_flex = report.plot_by_id["nexons_end_flex_plot"].datasets[0].lines[0]
    assert end_flex.pairs == [(x, x + 500) for x in range(-500, 501)]
