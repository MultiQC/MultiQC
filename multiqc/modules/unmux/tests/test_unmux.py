"""Tests for the unmux module.

Most tests run on `data/modules/unmux` in MultiQC/test-data: unmux 0.4.0 metrics from simulated reads for four
pools. Pool run1 has samples with sub-samples and removed reads. Pool run2.fq takes unmux's default names with
no sample sheet, and MultiQC trims it to run2. Pool run3 passes every read through and writes an empty
per-sample file, and run4 has ten samples of uneven depth. Edge cases use small inline files.
"""

import shutil
from pathlib import Path
from typing import Any, Dict, List

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.unmux import MultiqcModule
from multiqc.plots.table_object import Cell, InputRow
from multiqc.types import Anchor, ColumnKey, SampleGroup
from multiqc.utils import testing

PER_SAMPLE_HEADER = "pool\tsample\tsub_sample\treads\tfrac_of_pool\tfrac_bases_unextracted\n"
SUMMARY_HEADER = (
    "pool\ttotal_reads\treads_assigned\treads_pass_through\treads_unassigned\treads_removed\tfrac_assigned"
    "\tfrac_bases_unextracted\tmean_read_len_by_input\n"
)


@pytest.fixture
def test_data_dir() -> Path:
    return testing.data_dir() / "modules" / "unmux"


@pytest.fixture
def run_unmux(tmp_path):
    """Factory: copy the given test-data files and write the given inline files to a temp dir, then run unmux."""
    config.reset()

    def _run(files: Dict[str, str], copy_from: Any = None) -> MultiqcModule:
        if copy_from is not None:
            for path in Path(copy_from).glob("*.tsv"):
                shutil.copy(path, tmp_path / path.name)
        for filename, content in files.items():
            (tmp_path / filename).write_text(content)
        report.reset()
        config.preserve_module_raw_data = True
        report.analysis_files = [tmp_path]
        report.search_files(["unmux"])
        return MultiqcModule()

    yield _run
    report.reset()
    config.reset()


@pytest.fixture
def all_pools(run_unmux, test_data_dir):
    return lambda: run_unmux({}, copy_from=test_data_dir)


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


def test_every_pool_and_sample_is_found(all_pools):
    module = all_pools()

    assert set(saved(module, "multiqc_unmux_summary")) == {"run1", "run2", "run3", "run4"}
    assert len(saved(module, "multiqc_unmux_per_sample")) == 17


def test_files_are_found_by_header_not_name(run_unmux, test_data_dir):
    module = run_unmux(
        {
            "anything.txt": (test_data_dir / "run1.unmux.summary.tsv").read_text(),
            "else.tsv": (test_data_dir / "run1.unmux.per_sample.tsv").read_text(),
        }
    )

    assert set(saved(module, "multiqc_unmux_summary")) == {"run1"}
    assert len(saved(module, "multiqc_unmux_per_sample")) == 4


def test_samples_are_named_with_their_sub_sample(all_pools):
    module = all_pools()

    names = set(saved(module, "multiqc_unmux_per_sample"))
    assert {"tumor (lib01)", "normal (lib02)", "TAGGCTTC", "sampleA", "sampleJ"} <= names
    assert saved(module, "multiqc_unmux_per_sample")["tumor (lib01)"]["pool"] == "run1"


def test_outcomes_add_up_to_the_pool(all_pools):
    all_pools()

    assert bar_values("unmux_outcomes", 0, "run1") == {
        "Assigned": 2848,
        "Passed through": 0,
        "Removed": 471,
        "Unassigned": 681,
    }
    assert bar_values("unmux_outcomes", 0, "run3")["Passed through"] == 1196


def test_reads_per_sample_heads_each_pool_with_its_leftovers(all_pools):
    all_pools()

    # Horizontal bar plots keep their samples bottom to top.
    plot: Any = report.plot_by_id[Anchor("unmux_reads_per_sample")]
    samples = list(reversed(plot.datasets[0].samples))
    assert samples[:4] == ["run1: unassigned", "run1: removed", "tumor (lib01)", "tumor (lib02)"]
    assert samples[6] == "run2: unassigned"
    assert "run2: removed" not in samples
    assert not any(s.startswith("run3") for s in samples)
    assert bar_values("unmux_reads_per_sample", 1, "tumor (lib01)")["Assigned"] == pytest.approx(19.5)
    assert bar_values("unmux_reads_per_sample", 1, "sampleJ")["Assigned"] == pytest.approx(0.0167)
    assert bar_values("unmux_reads_per_sample", 1, "run1: removed")["Removed"] == pytest.approx(11.775)


def test_metrics_table_nests_samples_under_their_pool(all_pools):
    all_pools()

    rows = table_rows("run1")
    assert [str(row.sample) for row in rows][:3] == ["run1", "tumor (lib01)", "tumor (lib02)"]
    pool = rows[0]
    assert raw(pool, "total_reads") == 4000
    assert raw(pool, "reads") == 2848
    assert raw(pool, "frac_unassigned") == pytest.approx(0.17025)
    assert raw(pool, "frac_removed") == pytest.approx(0.11775)

    pass_through = table_rows("run3")
    assert len(pass_through) == 1
    assert raw(pass_through[0], "reads") == 1196
    assert raw(pass_through[0], "frac_of_pool") == pytest.approx(0.797333)


def test_general_stats_groups_samples_under_their_pool(all_pools):
    all_pools()

    rows = general_stats_rows()
    run4 = rows[SampleGroup("run4")]
    assert str(run4[0].sample) == "run4"
    assert run4[0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.938333, abs=1e-5)
    assert run4[0].data[ColumnKey("reads")] == 5630
    assert len(run4) == 11
    assert {str(row.sample) for row in rows[SampleGroup("run1")][1:]} == {
        "tumor (lib01)",
        "tumor (lib02)",
        "normal (lib01)",
        "normal (lib02)",
    }
    assert rows[SampleGroup("run3")][0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.797333, abs=1e-5)


def test_general_stats_grouping_can_be_turned_off(all_pools):
    setattr(config, "unmux_config", {"general_stats_group_by_pool": False})
    all_pools()

    rows = general_stats_rows()
    assert SampleGroup("run1") not in rows
    assert rows[SampleGroup("tumor (lib01)")][0].data[ColumnKey("frac_of_pool")] == pytest.approx(0.195)


def test_general_stats_hides_reads_by_default(all_pools):
    all_pools()

    headers = {k: v for h in report.general_stats_headers.values() for k, v in h.items()}
    assert not headers[ColumnKey("frac_of_pool")].get("hidden")
    assert headers[ColumnKey("reads")]["hidden"] is True


def test_concatenated_files_hold_several_pools(run_unmux, test_data_dir):
    module = run_unmux(
        {
            "all.summary.tsv": "".join(p.read_text() for p in sorted(test_data_dir.glob("*.summary.tsv"))),
            "all.per_sample.tsv": "".join(p.read_text() for p in sorted(test_data_dir.glob("*.per_sample.tsv"))),
        }
    )

    assert set(saved(module, "multiqc_unmux_summary")) == {"run1", "run2", "run3", "run4"}
    assert len(saved(module, "multiqc_unmux_per_sample")) == 17


def test_a_name_shared_by_two_pools_gets_the_pool_in_front(run_unmux):
    module = run_unmux(
        {
            "p1.tsv": PER_SAMPLE_HEADER + "run1\tS1\t\t10\t0.5\t0\nrun1\tS2\t\t10\t0.5\t0\n",
            "p2.tsv": PER_SAMPLE_HEADER + "run2\tS1\t\t10\t1.0\t0\nrun2\trun1\t\t0\t0.0\t0\n",
        }
    )

    assert set(saved(module, "multiqc_unmux_per_sample")) == {"run1: S1", "S2", "run2: S1", "run2: run1"}


def test_per_sample_without_a_summary_still_reports(run_unmux, test_data_dir):
    module = run_unmux({"run4.tsv": (test_data_dir / "run4.unmux.per_sample.tsv").read_text()})

    assert "multiqc_unmux_summary" not in (module.saved_raw_data or {})
    assert "unmux_outcomes" not in report.plot_by_id
    pool = table_rows("run4")[0]
    assert raw(pool, "reads") == 5630
    assert raw(pool, "frac_of_pool") == pytest.approx(0.938333, abs=1e-5)


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


def test_ignored_pools_and_samples_are_dropped(all_pools):
    config.sample_names_ignore = ["run3", "sampleI", "sampleJ"]
    module = all_pools()

    assert "run3" not in saved(module, "multiqc_unmux_summary")
    assert not {"sampleI", "sampleJ"} & set(saved(module, "multiqc_unmux_per_sample"))


@pytest.mark.parametrize(
    "content",
    [
        SUMMARY_HEADER + "poolZ\t1\t2\n",
        SUMMARY_HEADER + "poolZ\t\t0\t0\t0\t0\t0\t0\t\n",
        PER_SAMPLE_HEADER + "poolZ\tS1\t\tmany\t0.5\t0\n",
    ],
)
def test_malformed_files_fail_loudly(run_unmux, content):
    with pytest.raises(ValueError, match="unmux: bad.tsv"):
        run_unmux({"bad.tsv": content})
