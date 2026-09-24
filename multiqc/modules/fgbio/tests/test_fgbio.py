"""Inline-fixture unit tests for the fgbio module's ClipBam parser.

The count fixtures are hand-built rather than copied from a run so that every derived
percentage can be checked against arithmetic done by hand in the test. They keep the
identities fgbio guarantees: `Pair` is the sum of `ReadOne` and `ReadTwo`, `All` is the
sum of `Fragment` and `Pair`, and the per-reason clipped-base counts add up to
`bases_clipped_post`.
"""

import pytest

from multiqc import config, report
from multiqc.modules.fgbio import MultiqcModule
from multiqc.types import ColumnKey, SampleGroup

CLIPBAM_HEADER = (
    "read_type\treads\treads_unmapped\treads_clipped_pre\treads_clipped_post\treads_clipped_five_prime"
    "\treads_clipped_three_prime\treads_clipped_overlapping\treads_clipped_extending\tbases\tbases_clipped_pre"
    "\tbases_clipped_post\tbases_clipped_five_prime\tbases_clipped_three_prime\tbases_clipped_overlapping"
    "\tbases_clipped_extending\n"
)

# Paired-end data: every read is paired, so the Fragment row is all zeros and All equals Pair.
# All: 2000 reads, 1220 clipped; 178000 aligned bases remain, 25000 were clipped (22000 for overlap).
CLIPBAM_PAIRED_TSV = (
    CLIPBAM_HEADER + "Fragment\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    "ReadOne\t1000\t0\t100\t600\t0\t0\t550\t0\t90000\t2000\t12000\t0\t0\t10000\t0\n"
    "ReadTwo\t1000\t0\t50\t620\t0\t0\t600\t0\t88000\t1000\t13000\t0\t0\t12000\t0\n"
    "Pair\t2000\t0\t150\t1220\t0\t0\t1150\t0\t178000\t3000\t25000\t0\t0\t22000\t0\n"
    "All\t2000\t0\t150\t1220\t0\t0\t1150\t0\t178000\t3000\t25000\t0\t0\t22000\t0\n"
)

# Fragment-only data: the paired rows are all zeros and All equals Fragment.
# All: 500 reads, 200 clipped; 40000 aligned bases remain, 1300 were clipped (1000 at the 5' end).
CLIPBAM_FRAGMENT_TSV = (
    CLIPBAM_HEADER + "Fragment\t500\t0\t20\t200\t200\t0\t0\t0\t40000\t300\t1300\t1000\t0\t0\t0\n"
    "ReadOne\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    "ReadTwo\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    "Pair\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    "All\t500\t0\t20\t200\t200\t0\t0\t0\t40000\t300\t1300\t1000\t0\t0\t0\n"
)

# An empty input BAM: every count is zero, so no percentage can be computed.
CLIPBAM_EMPTY_TSV = CLIPBAM_HEADER + "".join(
    f"{read_type}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    for read_type in ("Fragment", "ReadOne", "ReadTwo", "Pair", "All")
)

GROUPREADSBYUMI_TSV = (
    "family_size\tcount\tfraction\tfraction_gt_or_eq_family_size\n"
    "1\t7092682\t0.56602\t1\n"
    "2\t2005532\t0.16005\t0.43398\n"
    "3\t1039240\t0.08294\t0.27393\n"
)

ERRORRATEBYREADPOSITION_TSV = (
    "read_number\tposition\tbases_total\terrors\terror_rate\ta_to_c_error_rate\ta_to_g_error_rate"
    "\ta_to_t_error_rate\tc_to_a_error_rate\tc_to_g_error_rate\tc_to_t_error_rate\n"
    "1\t1\t1000\t10\t0.01\t0.001\t0.002\t0.001\t0.002\t0.002\t0.002\n"
    "1\t2\t1000\t20\t0.02\t0.002\t0.004\t0.002\t0.004\t0.004\t0.004\n"
)


@pytest.fixture
def run_fgbio_module(tmp_path):
    """Factory: write each `filename: content` pair to a temp dir and run the fgbio module on it."""

    def _run(files: dict):
        for filename, content in files.items():
            (tmp_path / filename).write_text(content)
        report.reset()
        report.analysis_files = [tmp_path]
        report.search_files(["fgbio"])
        return MultiqcModule()

    return _run


def _clipbam_table_rows(sample: str):
    """Return the grouped rows for `sample` from the ClipBam table, or None if absent."""
    for plot in report.plot_by_id.values():
        for dataset in getattr(plot, "datasets", []):
            dt = getattr(dataset, "dt", None)
            if dt is None or dt.id != "fgbio-clipbam-table":
                continue
            for section in dt.section_by_id.values():
                if SampleGroup(sample) in section.rows_by_sgroup:
                    return section.rows_by_sgroup[SampleGroup(sample)]
    return None


def _general_stats_for(sample: str):
    """Return (headers, row data) for `sample` from whichever general stats section holds it."""
    for section_key, rows_by_group in report.general_stats_data.items():
        if SampleGroup(sample) in rows_by_group:
            rows = rows_by_group[SampleGroup(sample)]
            assert len(rows) == 1
            return report.general_stats_headers[section_key], rows[0].data
    raise AssertionError(f"{sample!r} not found in general stats")


def test_general_stats_use_pre_clip_bases_as_denominator(run_fgbio_module):
    """`bases` counts only the aligned bases left after clipping, so the base percentage
    must divide by `bases + bases_clipped_post` (203000 here), not by `bases` alone."""
    run_fgbio_module({"sampleA.txt": CLIPBAM_PAIRED_TSV})

    headers, data = _general_stats_for("sampleA")
    assert data[ColumnKey("pct_bases_clipped")] == pytest.approx(100.0 * 25000 / (178000 + 25000))
    assert data[ColumnKey("pct_bases_clipped")] != pytest.approx(100.0 * 25000 / 178000)
    assert data[ColumnKey("pct_reads_clipped")] == pytest.approx(100.0 * 1220 / 2000)

    assert not headers[ColumnKey("pct_bases_clipped")].get("hidden")
    assert headers[ColumnKey("pct_reads_clipped")]["hidden"] is True
    assert "bases + bases_clipped_post" in headers[ColumnKey("pct_bases_clipped")]["description"]


def test_table_groups_all_as_primary_and_drops_empty_fragment_row(run_fgbio_module):
    """The `All` row heads each sample's group under the bare sample name, with `Pair`,
    `ReadOne` and `ReadTwo` nested under it in that order. The all-zero `Fragment` row of
    paired-end data is left out."""
    run_fgbio_module({"sampleA.txt": CLIPBAM_PAIRED_TSV})

    rows = _clipbam_table_rows("sampleA")
    assert rows is not None, "ClipBam table did not produce a sampleA group"
    assert [str(r.sample) for r in rows] == ["sampleA", "sampleA (Pair)", "sampleA (ReadOne)", "sampleA (ReadTwo)"]

    assert rows[0].data[ColumnKey("reads")].raw == 2000
    assert rows[1].data[ColumnKey("reads")].raw == 2000
    assert rows[2].data[ColumnKey("reads")].raw == 1000
    assert rows[3].data[ColumnKey("reads")].raw == 1000
    assert rows[0].data[ColumnKey("bases_pre_clip")].raw == 178000 + 25000
    assert rows[2].data[ColumnKey("pct_bases_clipped_overlapping")].raw == pytest.approx(
        100.0 * 10000 / (90000 + 12000)
    )


def test_fragment_only_data_collapses_to_a_single_headline_row(run_fgbio_module):
    """With only unpaired reads `All` equals `Fragment`, so nesting the `Fragment` row would
    just repeat the headline; the group degrades to the headline row alone."""
    run_fgbio_module({"sampleB.txt": CLIPBAM_FRAGMENT_TSV})

    rows = _clipbam_table_rows("sampleB")
    assert rows is not None, "ClipBam table did not produce a sampleB group"
    assert [str(r.sample) for r in rows] == ["sampleB"]
    assert rows[0].data[ColumnKey("reads")].raw == 500
    assert rows[0].data[ColumnKey("pct_bases_clipped")].raw == pytest.approx(100.0 * 1300 / (40000 + 1300))

    _, data = _general_stats_for("sampleB")
    assert data[ColumnKey("pct_reads_clipped")] == pytest.approx(40.0)


def test_all_zero_metrics_add_no_rows_but_do_not_fail(run_fgbio_module):
    """An empty input BAM yields all-zero counts: nothing can be divided, so the sample gets
    no percentages, no table row and no general stats entry, without raising."""
    module = run_fgbio_module({"sampleC.txt": CLIPBAM_EMPTY_TSV})

    assert "fgbio-clipbam" not in {section.anchor for section in module.sections}
    assert _clipbam_table_rows("sampleC") is None
    assert not any(SampleGroup("sampleC") in rows for rows in report.general_stats_data.values())


def test_data_file_flattens_read_types_and_keeps_derived_metrics(run_fgbio_module):
    original = config.preserve_module_raw_data
    config.preserve_module_raw_data = True
    try:
        module = run_fgbio_module({"sampleA.txt": CLIPBAM_PAIRED_TSV})
    finally:
        config.preserve_module_raw_data = original

    saved = module.saved_raw_data
    assert saved is not None
    row = saved["multiqc_fgbio_clipbam"]["sampleA"]
    assert row["All_reads"] == 2000
    assert row["ReadTwo_bases_clipped_overlapping"] == 12000
    assert row["All_bases_pre_clip"] == 203000
    assert row["All_pct_bases_clipped"] == pytest.approx(100.0 * 25000 / 203000)
    # A zero denominator leaves the percentage out rather than reporting 0.
    assert "Fragment_pct_reads_clipped" not in row
    assert "Fragment_pct_bases_clipped" not in row


def test_clipbam_pattern_does_not_match_other_fgbio_outputs(run_fgbio_module):
    module = run_fgbio_module(
        {
            "sampleD.histo.tsv": GROUPREADSBYUMI_TSV,
            "sampleD.error_rate_by_read_position.txt": ERRORRATEBYREADPOSITION_TSV,
        }
    )

    assert list(module.find_log_files("fgbio/clipbam")) == []
    anchors = {section.anchor for section in module.sections}
    assert "fgbio-clipbam" not in anchors
    assert {"fgbio-groupreadsbyumi", "fgbio-error-rate-by-read-position"} <= anchors


def test_other_fgbio_patterns_do_not_match_clipbam_output(run_fgbio_module):
    module = run_fgbio_module({"sampleA.txt": CLIPBAM_PAIRED_TSV})

    assert list(module.find_log_files("fgbio/groupreadsbyumi")) == []
    assert list(module.find_log_files("fgbio/errorratebyreadposition")) == []
    assert {section.anchor for section in module.sections} == {"fgbio-clipbam-bases", "fgbio-clipbam"}
