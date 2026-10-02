import logging

import pytest

from multiqc import report
from multiqc.modules.mosdepth import MultiqcModule

GOOD_DIST = """\
chr1\t2\t0.90
chr1\t1\t0.95
chr1\t0\t1.00
total\t2\t0.90
total\t1\t0.95
total\t0\t1.00
"""

GOOD_SUMMARY = """\
chrom\tlength\tbases\tmean\tmin\tmax
chr1\t100\t500\t5.00\t0\t10
total\t100\t500\t5.00\t0\t10
"""


@pytest.mark.parametrize(
    "bad_fn,bad_content",
    [
        # Truncated mid-line, as reported in https://github.com/MultiQC/MultiQC/issues/3690
        ("bad.mosdepth.global.dist.txt", "chr1\t2\t0.90\nchr1\t1\t0.95\nchrUn_KI270438v1"),
        # Truncated cleanly before the "total" rows were written
        ("bad.mosdepth.global.dist.txt", "chr1\t2\t0.90\nchr1\t1\t0.95\nchr1\t0\t1.00\n"),
        # Non-numeric value
        ("bad.mosdepth.global.dist.txt", "chr1\t2\t0.90\ntotal\tx\t0.95\n"),
        # Truncated summary
        ("bad.mosdepth.summary.txt", "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t500"),
    ],
    ids=["dist_partial_line", "dist_no_total", "dist_non_numeric", "summary_partial_line"],
)
def test_malformed_file_is_skipped(tmp_path, caplog, bad_fn, bad_content):
    """A malformed file is skipped with a warning instead of crashing the module."""
    good_dist = tmp_path / "good.mosdepth.global.dist.txt"
    good_dist.write_text(GOOD_DIST)
    good_summary = tmp_path / "good.mosdepth.summary.txt"
    good_summary.write_text(GOOD_SUMMARY)
    bad_file = tmp_path / bad_fn
    bad_file.write_text(bad_content)

    report.reset()
    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])
    with caplog.at_level(logging.WARNING):
        MultiqcModule()

    assert f"Skipping {bad_file}" in caplog.text
    samples = {str(sample) for section in report.general_stats_data.values() for sample in section}
    assert samples == {"good"}
