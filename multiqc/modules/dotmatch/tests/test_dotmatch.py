from io import StringIO

from multiqc.modules.dotmatch.dotmatch import _parse_crispr_qc, _parse_sample_qc


def test_parse_sample_qc_preserves_ambiguity_and_rates():
    rows = _parse_sample_qc(
        StringIO(
            "sample_id\tassigned_reads\tambiguous_reads\tassignment_rate\tambiguous_rate\nsample_a\t8\t2\t0.8\t0.2\n"
        ),
        "fallback",
    )

    assert rows == {
        "sample_a": {
            "total_reads": "",
            "valid_extracted_reads": "",
            "assigned_reads": 8,
            "exact_reads": "",
            "k1_rescued_reads": "",
            "ambiguous_reads": 2,
            "no_match_reads": "",
            "invalid_reads": "",
            "assignment_rate": 0.8,
            "exact_rate": "",
            "rescue_rate": "",
            "ambiguous_rate": 0.2,
            "no_match_rate": "",
            "targets_observed": "",
            "zero_count_targets": "",
            "gini_index": "",
            "top_1pct_read_fraction": "",
            "candidates_verified": "",
        }
    }


def test_parse_crispr_qc_coerces_library_metrics():
    rows = _parse_crispr_qc(
        StringIO(
            "sample_id\tqc_status\ttotal_count\tcoverage_fraction\tzero_count_fraction\n"
            "sample_a\tpass\t12\t0.75\t0.25\n"
        ),
        "fallback",
    )

    assert rows["sample_a"]["qc_status"] == "pass"
    assert rows["sample_a"]["total_count"] == 12
    assert rows["sample_a"]["coverage_fraction"] == 0.75
    assert rows["sample_a"]["zero_count_fraction"] == 0.25
