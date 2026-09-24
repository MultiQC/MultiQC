import os
from pathlib import Path

import math
import pytest

from multiqc import config, report
from multiqc.modules.samtools import MultiqcModule
from multiqc.types import ColumnKey, SampleGroup, SectionKey
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


def test_data_parsed(data_dir):
    data_subdir = data_dir / "modules/samtools/flagstat"
    for path in os.listdir(data_subdir):
        path = data_subdir / path
        report.analysis_files = [path]
        report.search_files(["samtools"])
        config.preserve_module_raw_data = True
        m = MultiqcModule()
        assert m.saved_raw_data is not None
        assert len(m.saved_raw_data) > 0
        assert m._clean_s_name(Path(path).name) in list(m.saved_raw_data.values())[0]


def slurp_file(data_dir, fname):
    with (data_dir / "modules/samtools/flagstat" / fname).open() as fh:
        return fh.read()


@pytest.fixture
def rep1(data_dir):
    return slurp_file(data_dir, "small.samtools13.flagstat.log.txt")


@pytest.fixture
def rep2(data_dir):
    return slurp_file(data_dir, "small.samtools12.flagstat.log.txt")


def test_rep1(rep1):
    """Test that parsing rep1 produces expected results"""
    from multiqc.modules.samtools.flagstat import parse_single_report

    res1 = parse_single_report(rep1)

    # I expect 13 + 13 + 3 + 3 + 1 things reported in total
    assert len(res1) == 13 + 13 + 3 + 3 + 1

    assert (res1["total_passed"], res1["total_failed"]) == (5414, 0)

    assert res1["flagstat_total"] == 5414

    assert res1["mapped_passed_pct"] == 98.82

    # I expect mapped_failed_pct to be float('nan')
    assert math.isnan(res1["mapped_failed_pct"])


def test_rep2(rep1, rep2):
    """I expect rep2 to give identical results to rep1."""
    from multiqc.modules.samtools.flagstat import parse_single_report

    res1 = parse_single_report(rep1)
    res2 = parse_single_report(rep2)

    # But since nan != nan we have to strip these out.
    nans = [k for k, v in res1.items() if math.isnan(v)]
    for k in nans:
        del res1[k]
        del res2[k]

    assert res1 == res2


def _flagstat_report(total_passed: int, mapped_passed: int, total_failed: int, mapped_failed: int) -> str:
    mapped_passed_pct = mapped_passed / total_passed * 100
    mapped_failed_pct = mapped_failed / total_failed * 100
    return f"""\
{total_passed} + {total_failed} in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
{mapped_passed} + {mapped_failed} mapped ({mapped_passed_pct:.2f}% : {mapped_failed_pct:.2f}%)
{total_passed} + {total_failed} paired in sequencing
{total_passed // 2} + {total_failed // 2} read1
{total_passed // 2} + {total_failed // 2} read2
{mapped_passed} + {mapped_failed} properly paired ({mapped_passed_pct:.2f}% : {mapped_failed_pct:.2f}%)
{mapped_passed} + {mapped_failed} with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
"""


def test_grouped_lanes_sum_mapped_reads_in_general_stats(tmp_path):
    """Grouped lane-level flagstat reports should sum mapped reads."""
    report.reset()
    config.reset()

    (tmp_path / "sample_L001.flagstat").write_text(
        _flagstat_report(total_passed=100, mapped_passed=80, total_failed=20, mapped_failed=10)
    )
    (tmp_path / "sample_L002.flagstat").write_text(
        _flagstat_report(total_passed=200, mapped_passed=150, total_failed=100, mapped_failed=50)
    )

    config.table_sample_merge = {
        "L001": "_L001",
        "L002": "_L002",
    }
    report.analysis_files = [tmp_path]
    report.search_files(["samtools"])

    MultiqcModule()

    group_rows = report.general_stats_data[SectionKey("samtools")][SampleGroup("sample")]
    assert group_rows[0].data[ColumnKey("mapped_passed")] == 230
    assert group_rows[0].data[ColumnKey("total_passed")] == 300
    assert group_rows[0].data[ColumnKey("flagstat_total")] == 420
    assert group_rows[0].data[ColumnKey("mapped_passed_pct")] == pytest.approx(76.67, abs=0.01)
