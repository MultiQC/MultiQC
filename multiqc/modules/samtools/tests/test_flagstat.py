import os
from pathlib import Path

import math
import pytest

from multiqc import report
from multiqc.modules.samtools import MultiqcModule
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
        m = MultiqcModule()
        assert len(m.saved_raw_data) > 0
        assert m.clean_s_name(Path(path).name) in list(m.saved_raw_data.values())[0]


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
