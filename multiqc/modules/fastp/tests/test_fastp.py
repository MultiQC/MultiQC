import json

import pytest

from multiqc import report, reset
from multiqc.modules.fastp import MultiqcModule


@pytest.fixture(autouse=True)
def reset_multiqc_state():
    reset()
    yield
    reset()


def run_fastp_module(tmp_path) -> MultiqcModule:
    summary = {
        "total_reads": 100,
        "total_bases": 1000,
        "q30_bases": 900,
        "q30_rate": 0.9,
        "gc_content": 0.4,
    }
    (tmp_path / "sample.fastp.json").write_text(
        json.dumps(
            {
                "command": "fastp -i sample_R1.fastq.gz -o out.fastq.gz",
                "summary": {
                    "fastp_version": "0.23.4",
                    "before_filtering": {**summary, "total_reads": 100, "total_bases": 1000, "q30_bases": 900},
                    "after_filtering": {**summary, "total_reads": 80, "total_bases": 700, "q30_bases": 650},
                },
                "filtering_result": {"passed_filter_reads": 80},
                "duplication": {"rate": 0.0},
            }
        )
    )
    report.analysis_files = [tmp_path]
    report.search_files(["fastp"])
    return MultiqcModule()


def test_general_stats_reads_before_filtering(tmp_path):
    module = run_fastp_module(tmp_path)
    (sample,) = module.fastp_data.values()
    assert sample["before_filtering_total_reads"] == 100
    headers = {k for hdrs in report.general_stats_headers.values() for k in hdrs}
    assert "before_filtering_total_reads" in headers
    assert "filtering_result_passed_filter_reads" in headers
