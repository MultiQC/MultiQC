import shutil

import pytest

from multiqc import report
from multiqc.modules.fastqc.fastqc import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


def _run_fastqc(tmp_path, data_dir):
    shutil.copy(data_dir / "modules/fastqc/fastqc_data.txt", tmp_path / "fastqc_data.txt")
    report.reset()
    report.analysis_files = [tmp_path]
    report.search_files(["fastqc"])
    return MultiqcModule()


def _general_stats_value(key):
    rows = [row for groups in report.general_stats_data.values() for rows in groups.values() for row in rows]
    assert len(rows) == 1
    return rows[0].data[key]


def test_percent_fails_ignores_basic_statistics(tmp_path, data_dir):
    _run_fastqc(tmp_path, data_dir)
    # 11 checks other than Basic Statistics, of which only Kmer Content fails
    assert _general_stats_value("percent_fails") == pytest.approx(100 / 11)
