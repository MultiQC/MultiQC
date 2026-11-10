"""Inline-fixture unit tests for the unmux module.

The fixtures are unmux 0.4.0 output for synthetic pools: one split into samples with two sub-samples each and a
`--remove` rule, one split into samples with no sub-samples, and one run with no samples at all, which passes
every read through and writes an empty per-sample file.
"""

from typing import Any, Dict, List

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.unmux import MultiqcModule
from multiqc.plots.table_object import Cell, InputRow
from multiqc.types import Anchor, ColumnKey, SampleGroup

PER_SAMPLE_HEADER = "pool\tsample\tsub_sample\treads\tfrac_of_pool\tfrac_bases_unextracted\n"
SUMMARY_HEADER = (
    "pool\ttotal_reads\treads_assigned\treads_pass_through\treads_unassigned\treads_removed\tfrac_assigned"
    "\tfrac_bases_unextracted\tmean_read_len_by_input\n"
)

POOL_A_PER_SAMPLE = (
    PER_SAMPLE_HEADER + "poolA\ttumor\tlib1\t994\t0.198800\t0.189189\n"
    "poolA\ttumor\tlib2\t1032\t0.206400\t0.189189\n"
    "poolA\tnormal\tlib1\t609\t0.121800\t0.189189\n"
    "poolA\tnormal\tlib2\t624\t0.124800\t0.189189\n"
    "poolA\tblank\tlib1\t194\t0.038800\t0.189189\n"
    "poolA\tblank\tlib2\t199\t0.039800\t0.189189\n"
)
POOL_A_SUMMARY = SUMMARY_HEADER + "poolA\t5000\t3652\t0\t848\t500\t0.730400\t0.189189\t74\n"

POOL_B_PER_SAMPLE = (
    PER_SAMPLE_HEADER + "poolB\tS1\t\t892\t0.297333\t0.117647\n"
    "poolB\tS2\t\t1807\t0.602333\t0.117647\n"
    "poolB\tS3\t\t7\t0.002333\t0.117647\n"
    "poolB\tCTL\t\t0\t0.000000\t0.000000\n"
)
POOL_B_SUMMARY = SUMMARY_HEADER + "poolB\t3000\t2706\t0\t294\t0\t0.902000\t0.117647\t68\n"

POOL_C_SUMMARY = SUMMARY_HEADER + "poolC\t2000\t0\t1696\t0\t304\t0.000000\t0.117647\t68\n"


@pytest.fixture
def run_unmux(tmp_path):
    """Factory: write each `filename: content` pair to a temp dir and run the unmux module on it."""
    config.reset()

    def _run(files: Dict[str, str]) -> MultiqcModule:
        for filename, content in files.items():
            (tmp_path / filename).write_text(content)
        report.reset()
        config.preserve_module_raw_data = True
        report.analysis_files = [tmp_path]
        report.search_files(["unmux"])
        return MultiqcModule()

    yield _run
    config.reset()


def all_pools(run_unmux) -> MultiqcModule:
    return run_unmux(
        {
            "poolA.unmux.per_sample.tsv": POOL_A_PER_SAMPLE,
            "poolA.unmux.summary.tsv": POOL_A_SUMMARY,
            "poolB.unmux.per_sample.tsv": POOL_B_PER_SAMPLE,
            "poolB.unmux.summary.tsv": POOL_B_SUMMARY,
            "poolC.unmux.per_sample.tsv": "",
            "poolC.unmux.summary.tsv": POOL_C_SUMMARY,
        }
    )


def general_stats_rows() -> Dict[SampleGroup, List[InputRow]]:
    rows: Dict[SampleGroup, List[InputRow]] = {}
    for rows_by_group in report.general_stats_data.values():
        rows.update(rows_by_group)
    return rows


def table_rows(group: str) -> List[InputRow]:
    plot: Any = report.plot_by_id[Anchor("unmux_metrics_table")]
    for dataset in plot.datasets:
        for section in dataset.dt.section_by_id.values():
            if SampleGroup(group) in section.rows_by_sgroup:
                return section.rows_by_sgroup[SampleGroup(group)]
    raise AssertionError(f"{group!r} not in the metrics table")


def raw(row: InputRow, column: str) -> Any:
    cell = row.data[ColumnKey(column)]
    return cell.raw if isinstance(cell, Cell) else cell


def saved(module: MultiqcModule, key: str) -> Dict[str, Any]:
    assert module.saved_raw_data is not None
    return module.saved_raw_data[key]


def bar_values(plot_id: str, dataset_idx: int, sample: str) -> Dict[str, float]:
    plot: Any = report.plot_by_id[Anchor(plot_id)]
    dataset = plot.datasets[dataset_idx]
    index = dataset.samples.index(sample)
    return {cat.name: cat.data[index] for cat in dataset.cats if cat.data[index] == cat.data[index]}


def test_no_files_raises(run_unmux):
    with pytest.raises(ModuleNoSamplesFound):
        run_unmux({"other.tsv": "a\tb\n1\t2\n"})


def test_files_are_found_by_header_not_name(run_unmux):
    module = run_unmux({"anything.txt": POOL_A_SUMMARY, "else.tsv": POOL_A_PER_SAMPLE})

    assert set(saved(module, "multiqc_unmux_summary")) == {"poolA"}
    assert len(saved(module, "multiqc_unmux_per_sample")) == 6


def test_samples_are_named_with_their_sub_sample(run_unmux):
    module = all_pools(run_unmux)

    names = set(saved(module, "multiqc_unmux_per_sample"))
    assert {"tumor (lib1)", "tumor (lib2)", "blank (lib2)", "S1", "CTL"} <= names
    assert saved(module, "multiqc_unmux_per_sample")["tumor (lib1)"]["pool"] == "poolA"


def test_outcomes_add_up_to_the_pool(run_unmux):
    all_pools(run_unmux)

    assert bar_values("unmux_outcomes", 0, "poolA") == {
        "Assigned": 3652,
        "Passed through": 0,
        "Removed": 500,
        "Unassigned": 848,
    }
    assert bar_values("unmux_outcomes", 0, "poolC")["Passed through"] == 1696


def test_reads_per_sample_heads_each_pool_with_its_leftovers(run_unmux):
    all_pools(run_unmux)

    # Horizontal bar plots keep their samples bottom to top.
    plot: Any = report.plot_by_id[Anchor("unmux_reads_per_sample")]
    samples = list(reversed(plot.datasets[0].samples))
    assert samples[:4] == ["poolA: unassigned", "poolA: removed", "tumor (lib1)", "tumor (lib2)"]
    assert samples[8] == "poolB: unassigned"
    assert "poolB: removed" not in samples
    assert not any(s.startswith("poolC") for s in samples)
    assert bar_values("unmux_reads_per_sample", 1, "tumor (lib1)")["Assigned"] == pytest.approx(19.88)
    assert bar_values("unmux_reads_per_sample", 1, "poolA: removed")["Removed"] == pytest.approx(10.0)


def test_metrics_table_nests_samples_under_their_pool(run_unmux):
    all_pools(run_unmux)

    rows = table_rows("poolA")
    assert [str(row.sample) for row in rows][:3] == ["poolA", "tumor (lib1)", "tumor (lib2)"]
    pool = rows[0]
    assert raw(pool, "total_reads") == 5000
    assert raw(pool, "reads") == 3652
    assert raw(pool, "frac_unassigned") == pytest.approx(0.1696)
    assert raw(pool, "frac_removed") == pytest.approx(0.1)

    pass_through = table_rows("poolC")
    assert len(pass_through) == 1
    assert raw(pass_through[0], "reads") == 1696
    assert raw(pass_through[0], "frac_of_pool") == pytest.approx(0.848)


def test_general_stats_groups_samples_under_their_pool(run_unmux):
    all_pools(run_unmux)

    rows = general_stats_rows()
    pool_a = rows[SampleGroup("poolA")]
    assert str(pool_a[0].sample) == "poolA"
    assert pool_a[0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.7304)
    assert pool_a[0].data[ColumnKey("reads")] == 3652
    assert {str(row.sample) for row in pool_a[1:]} == {
        "tumor (lib1)",
        "tumor (lib2)",
        "normal (lib1)",
        "normal (lib2)",
        "blank (lib1)",
        "blank (lib2)",
    }
    assert rows[SampleGroup("poolC")][0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.848)


def test_general_stats_grouping_can_be_turned_off(run_unmux):
    setattr(config, "unmux_config", {"general_stats_group_by_pool": False})
    all_pools(run_unmux)

    rows = general_stats_rows()
    assert SampleGroup("poolA") not in rows
    assert rows[SampleGroup("tumor (lib1)")][0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.1988)


def test_general_stats_hides_reads_by_default(run_unmux):
    all_pools(run_unmux)

    headers = {k: v for h in report.general_stats_headers.values() for k, v in h.items()}
    assert not headers[ColumnKey("frac_of_pool")].get("hidden")
    assert headers[ColumnKey("reads")]["hidden"] is True


def test_concatenated_files_hold_several_pools(run_unmux):
    module = run_unmux(
        {
            "all.summary.tsv": POOL_A_SUMMARY + POOL_B_SUMMARY + POOL_C_SUMMARY,
            "all.per_sample.tsv": POOL_A_PER_SAMPLE + POOL_B_PER_SAMPLE,
        }
    )

    assert set(saved(module, "multiqc_unmux_summary")) == {"poolA", "poolB", "poolC"}
    assert len(saved(module, "multiqc_unmux_per_sample")) == 10


def test_a_name_shared_by_two_pools_gets_the_pool_in_front(run_unmux):
    module = run_unmux(
        {
            "p1.tsv": PER_SAMPLE_HEADER + "run1\tS1\t\t10\t0.5\t0\nrun1\tS2\t\t10\t0.5\t0\n",
            "p2.tsv": PER_SAMPLE_HEADER + "run2\tS1\t\t10\t1.0\t0\nrun2\trun1\t\t0\t0.0\t0\n",
        }
    )

    assert set(saved(module, "multiqc_unmux_per_sample")) == {"run1: S1", "S2", "run2: S1", "run2: run1"}


def test_per_sample_without_a_summary_still_reports(run_unmux):
    module = run_unmux({"poolB.unmux.per_sample.tsv": POOL_B_PER_SAMPLE})

    assert "multiqc_unmux_summary" not in (module.saved_raw_data or {})
    assert "unmux_outcomes" not in report.plot_by_id
    pool = table_rows("poolB")[0]
    assert raw(pool, "reads") == 2706
    assert raw(pool, "frac_of_pool") == pytest.approx(0.901999)


def test_extra_columns_are_kept(run_unmux):
    per_sample = (
        "pool\tsample\tsub_sample\treads\tfrac_of_pool\tfrac_bases_unextracted\tmean_errors\tbarcode\n"
        "poolX\tS1\t\t10\t0.5\t0.1\t0.25\tACGT\n"
        "poolX\tS2\t\t10\t0.5\t0.1\t0.5\tTTTT\n"
    )
    module = run_unmux({"x.tsv": per_sample})

    assert saved(module, "multiqc_unmux_per_sample")["S1"]["mean_errors"] == 0.25
    sample = table_rows("poolX")[1]
    assert raw(sample, "mean_errors") == 0.25
    assert raw(sample, "barcode") == "ACGT"


def test_ignored_pools_and_samples_are_dropped(run_unmux):
    config.sample_names_ignore = ["poolC", "blank*"]
    module = all_pools(run_unmux)

    assert "poolC" not in saved(module, "multiqc_unmux_summary")
    assert not any(s.startswith("blank") for s in saved(module, "multiqc_unmux_per_sample"))


def test_malformed_rows_are_skipped(run_unmux):
    module = run_unmux({"a.tsv": POOL_A_SUMMARY + "poolZ\t1\t2\n" + "poolY\t\t0\t0\t0\t0\t0\t0\t\n"})

    assert set(saved(module, "multiqc_unmux_summary")) == {"poolA"}
