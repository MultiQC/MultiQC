import json
import logging
from typing import Any, Dict, List, Optional

import pytest

from multiqc import report, reset
from multiqc.modules.fastp import MultiqcModule


@pytest.fixture(autouse=True)
def reset_multiqc_state():
    reset()
    yield
    reset()


def run_fastp_module(tmp_path, histogram: Optional[List[int]], unknown: Optional[int]) -> MultiqcModule:
    insert_size: Dict[str, Any] = {"peak": 2}
    if histogram is not None:
        insert_size["histogram"] = histogram
    if unknown is not None:
        insert_size["unknown"] = unknown

    report_path = tmp_path / "sample.fastp.json"
    report_path.write_text(
        json.dumps(
            {
                "summary": {
                    "fastp_version": "0.23.4",
                    "before_filtering": {
                        "total_reads": 100,
                        "total_bases": 1000,
                    },
                    "after_filtering": {
                        "total_reads": 100,
                        "total_bases": 1000,
                    },
                },
                "filtering_result": {
                    "passed_filter_reads": 100,
                },
                "duplication": {
                    "rate": 0.0,
                    "histogram": [100],
                },
                "insert_size": insert_size,
                "command": "fastp --in1 sample.fastq.gz",
            }
        )
    )

    report.analysis_files = [report_path]
    report.search_files(["fastp"])
    return MultiqcModule()


@pytest.mark.parametrize(
    ("histogram", "unknown", "expected"),
    [
        ([10, 20, 0], 70, {1: 10.0, 2: 20.0}),
        ([10, 20, 0], 0, {1: 100.0 / 3.0, 2: 200.0 / 3.0}),
        ([10, 30, 0], 0, {1: 25.0, 2: 75.0}),
    ],
)
def test_insert_size_percentages_include_unknown_pairs(
    tmp_path,
    histogram: List[int],
    unknown: int,
    expected: Dict[int, float],
):
    module = run_fastp_module(tmp_path, histogram, unknown)

    assert module.fastp_insert_size_data["sample"] == pytest.approx(expected)


@pytest.mark.parametrize("unknown", [70, 0])
def test_insert_size_omitted_without_classified_pairs(tmp_path, unknown: int):
    module = run_fastp_module(tmp_path, [0, 0, 0], unknown)

    assert module.fastp_insert_size_data == {}


@pytest.mark.parametrize(
    ("histogram", "unknown"),
    [
        (None, 70),
        ([10, 20, 0], None),
    ],
)
def test_insert_size_omitted_when_required_count_is_missing(tmp_path, caplog, histogram, unknown):
    with caplog.at_level(logging.DEBUG, logger="multiqc.modules.fastp.fastp"):
        module = run_fastp_module(tmp_path, histogram, unknown)

    assert module.fastp_insert_size_data == {}
    assert "No insert size plot data: sample" in caplog.text
