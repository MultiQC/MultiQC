"""Tests for the fgmetric-shaped schema base class."""

import math
from typing import Optional

import pytest

from multiqc.modules.fgumi.schemas import FamilySizeMetric, FgumiMetric, MetricFormatError


class _Probe(FgumiMetric):
    name: str
    count: int
    rate: float
    size: Optional[int]


def test_read_parses_rows_in_order():
    rows = _Probe.read("name\tcount\trate\tsize\na\t1\t0.5\t7\nb\t2\t1\t\n", "probe.txt")
    assert [(r.name, r.count, r.rate, r.size) for r in rows] == [("a", 1, 0.5, 7), ("b", 2, 1.0, None)]


@pytest.mark.parametrize(
    "token,check",
    [("NaN", math.isnan), ("Infinity", lambda x: x == math.inf), ("-Infinity", lambda x: x == -math.inf)],
)
def test_read_accepts_fgbio_non_finite_tokens(token, check):
    (row,) = _Probe.read(f"name\tcount\trate\tsize\na\t1\t{token}\t\n", "probe.txt")
    assert check(row.rate)


def test_read_ignores_extra_columns_and_blank_lines():
    (row,) = _Probe.read("extra\tname\tcount\trate\tsize\n\nx\ta\t1\t0.5\t2\n\n", "probe.txt")
    assert row.name == "a"


def test_missing_column_names_file():
    with pytest.raises(MetricFormatError, match=r"probe\.txt.*missing.*rate"):
        _Probe.read("name\tcount\tsize\na\t1\t2\n", "probe.txt")


def test_ragged_row_names_file_and_line():
    with pytest.raises(MetricFormatError, match=r"probe\.txt: line 3"):
        _Probe.read("name\tcount\trate\tsize\na\t1\t0.5\t2\nb\t2\n", "probe.txt")


def test_bad_value_names_file():
    with pytest.raises(MetricFormatError, match=r"probe\.txt: data row 1"):
        _Probe.read("name\tcount\trate\tsize\na\tNOT_AN_INT\t0.5\t2\n", "probe.txt")


def test_empty_file_is_an_error():
    with pytest.raises(MetricFormatError, match="empty"):
        _Probe.read("", "probe.txt")


def test_schema_columns_follow_declaration_order():
    assert FamilySizeMetric.columns() == ["family_size", "count", "fraction", "fraction_gt_or_eq_family_size"]
