from pathlib import Path

import pytest

from multiqc import config, report, validation
from multiqc.modules.fgbio import MultiqcModule


@pytest.fixture(autouse=True)
def reset():
    report.reset()
    config.reset()
    validation.reset()


@pytest.fixture
def data_dir():
    return Path(__file__).parent


def test_family_sizes_parsed(data_dir):
    report.analysis_files = [str(data_dir)]
    report.search_files(["fgbio"])
    m = MultiqcModule()
    assert len(m.sections) > 0


def test_family_sizes_general_stats(data_dir):
    report.analysis_files = [str(data_dir)]
    report.search_files(["fgbio"])
    MultiqcModule()
    assert "fgbio" in report.general_stats_data
    fgbio_stats = report.general_stats_data["fgbio"]
    assert len(fgbio_stats) == 2
    for s_name, rows in fgbio_stats.items():
        for row in rows:
            assert "duplex_rate" in row.data
            assert 0 <= row.data["duplex_rate"] <= 1
            assert "mean_family_size" in row.data
            assert row.data["mean_family_size"] > 0


def test_family_sizes_two_samples(data_dir):
    report.analysis_files = [str(data_dir)]
    report.search_files(["fgbio"])
    m = MultiqcModule()
    duplex_sections = [
        s for s in m.sections
        if "Family Sizes" in getattr(s, "name", "")
    ]
    assert len(duplex_sections) == 1
