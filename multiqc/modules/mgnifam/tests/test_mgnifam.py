import json

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.mgnifam import MultiqcModule
from multiqc.types import ColumnKey, SampleGroup, SectionKey


def _stats(command="generate_families", **changes):
    stats = {
        "tool": "mgnifam",
        "schema_version": 1,
        "version": "3.1.0",
        "command": command,
        "chunk_id": "1",
        "exit_status": 0,
        "parameters": {"discard_min_rep_length": 75},
        "families": {"input": 4, "successful": 3, "discarded": 1, "converged": 2, "crashed": 0},
        "discard_reasons": {"few seed sequences remained": 1},
        "histograms": {
            "full_msa_size": {"4": 2, "16": 1},
            "model_length": {"120": 3},
            "representative_length": {"117": 3},
        },
    }
    if command == "update_families":
        stats["histograms"].update(
            {"model_length_change": {"-1": 1, "0": 2}, "rounds_run": {"2": 3}, "retention": {"0.5": 1, "1.0": 2}}
        )
    stats.update(changes)
    return stats


@pytest.fixture
def run_module(tmp_path):
    def _run(files):
        for name, stats in files.items():
            (tmp_path / name).write_text(json.dumps(stats, indent=2) + "\n")
        report.reset()
        config.reset()
        report.analysis_files = [tmp_path]
        report.search_files(["mgnifam"])
        return MultiqcModule()

    return _run


def _general_stats(sample):
    return report.general_stats_data[SectionKey("mgnifam")][SampleGroup(sample)][0].data


def test_both_commands_parse_as_separate_samples(run_module):
    module = run_module({"1_stats.json": _stats(), "1_updated_stats.json": _stats("update_families")})

    assert set(module.mgnifam_data) == {"1", "1_updated"}
    generated = _general_stats("1")
    assert generated[ColumnKey("successful_pct")] == pytest.approx(75)
    assert ColumnKey("mean_retention") not in generated
    assert _general_stats("1_updated")[ColumnKey("mean_retention")] == pytest.approx(2.5 / 3)
    # Update-only sections are added once, alongside the shared ones.
    anchors = [section.anchor for section in module.sections]
    assert anchors == [
        "mgnifam-outcomes",
        "mgnifam-full-msa-size",
        "mgnifam-model-length",
        "mgnifam-model-length-change",
        "mgnifam-retention",
    ]


def test_crashed_chunks_are_flagged(run_module):
    crashed = _stats(
        exit_status=3, families={"input": 4, "successful": 3, "discarded": 1, "converged": 2, "crashed": 1}
    )
    module = run_module({"1_stats.json": crashed, "2_stats.json": _stats()})

    section = next(s for s in module.sections if s.anchor == "mgnifam-outcomes")
    assert section.alerts and section.alerts[0].affected_samples == ["1"]


def test_update_sections_are_skipped_without_update_samples(run_module):
    module = run_module({"1_stats.json": _stats()})
    anchors = [section.anchor for section in module.sections]
    assert "mgnifam-retention" not in anchors


def test_empty_chunk_renders_without_a_percentage(run_module):
    empty = _stats(
        families={"input": 0, "successful": 0, "discarded": 0, "converged": 0, "crashed": 0},
        discard_reasons={},
        histograms={"full_msa_size": {}, "model_length": {}, "representative_length": {}},
    )
    module = run_module({"1_stats.json": empty})

    assert ColumnKey("successful_pct") not in _general_stats("1")
    # The section stays, with the chunk listed as having nothing to plot.
    section = next(s for s in module.sections if s.anchor == "mgnifam-full-msa-size")
    assert section.plot_anchor is None


def test_unsupported_schema_is_skipped(run_module):
    with pytest.raises(ModuleNoSamplesFound):
        run_module({"1_stats.json": _stats(schema_version=2)})


def test_other_tools_stats_files_are_not_matched(run_module):
    with pytest.raises(ModuleNoSamplesFound):
        run_module({"1_stats.json": {"tool": "something_else", "families": {}}})


def test_invalid_and_schemaless_files_are_skipped(run_module, tmp_path):
    (tmp_path / "broken_stats.json").write_text('{\n  "tool": "mgnifam",\n  "schema_')
    (tmp_path / "old_stats.json").write_text('{\n  "tool": "mgnifam"\n}\n')
    module = run_module({"1_stats.json": _stats()})
    assert set(module.mgnifam_data) == {"1"}


def test_retention_is_binned_for_display(run_module):
    stats = _stats("update_families")
    stats["histograms"]["retention"] = {"0.91": 1, "0.93": 2, "1.0": 3}
    run_module({"1_updated_stats.json": stats})

    plot = next(p for p in report.plot_by_id.values() if not isinstance(p, str) and p.id == "mgnifam-retention-plot")
    points = dict(plot.datasets[0].lines[0].pairs)
    assert points == {0.9: 3, 1.0: 3}
