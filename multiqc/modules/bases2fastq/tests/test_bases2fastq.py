"""Tests for bases2fastq module: parsers and integration."""

import json
from pathlib import Path
from typing import Any, List
from unittest.mock import patch

import pytest

from multiqc import report, config
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.types import LoadedFileDict

from multiqc.modules.bases2fastq.bases2fastq import MultiqcModule, _get_min_polonies


def _load_fixture(fixtures_dir: Path, *parts: str) -> dict:
    """Load JSON fixture; path is fixtures_dir / path0 / path1 / ... / filename."""
    path = fixtures_dir.joinpath(*parts)
    with path.open() as f:
        return json.load(f)


class TestExtractManifestLaneSettings:
    """Tests for _extract_manifest_lane_settings helper."""

    def test_extract_manifest_lane_settings_minimal(self, fixtures_dir):
        """Manifest with one lane yields run_lane -> settings."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_manifest = _load_fixture(fixtures_dir, "PairedEndDefaultProject", "RunManifest.json")
        result = m._extract_manifest_lane_settings(run_manifest, "RUN01-a1b2")
        assert len(result) == 1
        run_lane = "RUN01-a1b2 | L1"
        assert run_lane in result
        assert result[run_lane]["AdapterTrimType"] == "Paired-End"
        assert result[run_lane]["R1AdapterMinimumTrimmedLength"] == 16
        assert result[run_lane]["R2AdapterMinimumTrimmedLength"] == 16
        assert "Indexing" in result[run_lane]

    def test_extract_manifest_lane_settings_empty_settings(self, fixtures_dir, tmp_path):
        """Manifest without Settings returns empty dict."""
        report.reset()
        run_stats = _load_fixture(fixtures_dir, "PairedEndNoProject", "RunStats.json")
        (tmp_path / "RunStats.json").write_text(json.dumps(run_stats))
        (tmp_path / "RunManifest.json").write_text(json.dumps({}))
        report.analysis_files = [str(tmp_path)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        result = m._extract_manifest_lane_settings({}, "RUN01-a1b2")
        assert result == {}


class TestBuildIndexAssignmentFromStats:
    """Tests for _build_index_assignment_from_stats helper."""

    def test_build_index_assignment_from_stats_with_occurrences(self, fixtures_dir):
        """Project RunStats with Occurrences produces run_inner and percentages."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        project_stats = _load_fixture(
            fixtures_dir, "PairedEndDefaultProject", "Samples", "DefaultProject", "DefaultProject_RunStats.json"
        )
        run_inner, total = m._build_index_assignment_from_stats(project_stats, "RUN01-a1b2", project="DefaultProject")
        assert total == 100000
        assert "AAATTT" in run_inner
        assert run_inner["AAATTT"]["SamplePolonyCounts"] == 5000
        assert run_inner["AAATTT"]["PercentOfPolonies"] == 5.0

    def test_build_index_assignment_from_stats_project(self, fixtures_dir):
        """Project-level stats (Samples/DefaultProject) add Project key to entries."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        project_stats = _load_fixture(
            fixtures_dir, "PairedEndDefaultProject", "Samples", "DefaultProject", "DefaultProject_RunStats.json"
        )
        run_inner, _ = m._build_index_assignment_from_stats(project_stats, "RUN01-a1b2", project="DefaultProject")
        assert run_inner
        for entry in run_inner.values():
            assert entry.get("Project") == "DefaultProject"


class TestMergeManifestIndexSequences:
    """Tests for _merge_manifest_index_sequences helper."""

    def test_merge_manifest_index_sequences(self, fixtures_dir):
        """Index1/Index2 from RunManifest Samples merged into assignment dict."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        project_stats = _load_fixture(
            fixtures_dir, "PairedEndDefaultProject", "Samples", "DefaultProject", "DefaultProject_RunStats.json"
        )
        run_manifest = _load_fixture(fixtures_dir, "PairedEndDefaultProject", "RunManifest.json")
        run_inner, _ = m._build_index_assignment_from_stats(project_stats, "RUN01-a1b2", project="DefaultProject")
        sample_to_index = {"RUN01-a1b2": run_inner}
        m._merge_manifest_index_sequences(sample_to_index, run_manifest, "RUN01-a1b2")
        assert run_inner["AAATTT"]["Index1"] == "AAA"
        assert run_inner["AAATTT"]["Index2"] == "TTT"


class TestParseRunProjectData:
    """Tests for run-level and project-level parsing."""

    def test_parse_run_project_data_run_level(self, fixtures_dir):
        """Run-level only (PairedEndNoProject) populates run and sample dicts."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert len(m.run_level_data) >= 1
        assert len(m.run_level_samples) >= 1
        sample_id = next(iter(m.run_level_samples))
        assert "__" in sample_id
        assert sample_id in m.run_level_samples_to_project

    def test_parse_run_project_data_min_polonies_filter(self, fixtures_dir):
        """Samples below min_polonies excluded (PairedEndNoProjectLowPolonies, config lowered)."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProjectLowPolonies"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        import multiqc.modules.bases2fastq.bases2fastq as b2f_mod

        with patch.object(b2f_mod, "_get_min_polonies", return_value=100):
            m = MultiqcModule()
        assert len(m.run_level_samples) == 1
        sample_id = next(iter(m.run_level_samples))
        assert sample_id.endswith("__Sample2")
        assert not any(s.endswith("__Sample1") for s in m.run_level_samples)


class TestParseRunUnassignedSequences:
    """Tests for unassigned sequences parser."""

    def test_parse_run_unassigned_sequences(self, fixtures_dir):
        """RunStats with Lanes/UnassignedSequences (PairedEndNoProjectWithLanes) produces int-keyed dict."""
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProjectWithLanes"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        unassigned = m._parse_run_unassigned_sequences("bases2fastq/run")
        assert isinstance(unassigned, dict)
        for k, v in unassigned.items():
            assert isinstance(k, int)
            assert "Run Name" in v
            assert "I1" in v
            assert "I2" in v
            assert "Number of Polonies" in v


class TestGetMinPolonies:
    """Tests for _get_min_polonies config helper."""

    def test_get_min_polonies_default_when_config_not_dict(self):
        with patch.object(config, "bases2fastq_config", None, create=True):
            assert _get_min_polonies() == 1000
        with patch.object(config, "bases2fastq_config", "string", create=True):
            assert _get_min_polonies() == 1000

    def test_get_min_polonies_invalid_int_uses_default(self):
        with patch.object(config, "bases2fastq_config", {"min_polonies": "bad"}, create=True):
            assert _get_min_polonies() == 1000
        with patch.object(config, "bases2fastq_config", {"min_polonies": None}, create=True):
            assert _get_min_polonies() == 1000

    def test_get_min_polonies_custom_value(self):
        with patch.object(config, "bases2fastq_config", {"min_polonies": 5000}, create=True):
            assert _get_min_polonies() == 5000


class TestValidatePath:
    """Tests for _validate_path security check."""

    def test_validate_path_escaped_returns_false(self, fixtures_dir, tmp_path):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        base = tmp_path / "sub"
        base.mkdir()
        outside = base.parent.parent.resolve()
        assert m._validate_path(outside / "any", base.resolve()) is False

    def test_validate_path_inside_returns_true(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._validate_path(run_dir / "RunStats.json", run_dir) is True


class TestReadJsonFile:
    """Tests for _read_json_file with validation and errors."""

    def test_read_json_file_path_outside_base_returns_none(self, fixtures_dir, tmp_path):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        outside = (tmp_path / "..").resolve()
        assert m._read_json_file(outside / "any.json", base_directory=tmp_path) is None

    def test_read_json_file_missing_file_returns_none(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._read_json_file(run_dir / "DoesNotExist.json") is None

    def test_read_json_file_invalid_json_returns_none(self, fixtures_dir, tmp_path):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        bad = tmp_path / "bad.json"
        bad.write_text("not json {")
        assert m._read_json_file(bad) is None


class TestExtractRunAnalysisName:
    """Tests for _extract_run_analysis_name."""

    def test_extract_run_analysis_name_missing_runname_returns_none(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._extract_run_analysis_name({"AnalysisID": "a1b2"}, "test") is None

    def test_extract_run_analysis_name_missing_analysisid_returns_none(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._extract_run_analysis_name({"RunName": "RUN01"}, "test") is None

    def test_extract_run_analysis_name_ok(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._extract_run_analysis_name({"RunName": "RUN01", "AnalysisID": "a1b2c3d4"}) == "RUN01-a1b2"


class TestParseRunProjectDataEdgeCases:
    """Edge cases for _parse_run_project_data."""

    def test_parse_run_project_data_empty_data_source_returns_empty(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_data, sample_data, sample_to_project = m._parse_run_project_data("", log_files=[])
        assert run_data == {}
        assert sample_data == {}
        assert sample_to_project == {}

    def test_parse_run_project_data_ignore_sample_skips_run(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        run_stats = _load_fixture(fixtures_dir, "PairedEndNoProject", "RunStats.json")
        log_files: List[LoadedFileDict[Any]] = [
            {
                "f": json.dumps(run_stats),
                "root": str(run_dir),
                "fn": "RunStats.json",
                "sp_key": "bases2fastq/run",
                "s_name": "RUN01-a1b2",
            }
        ]
        m = MultiqcModule()
        with patch.object(m, "is_ignore_sample", return_value=True):
            run_data, sample_data, _ = m._parse_run_project_data("bases2fastq/run", log_files=log_files)
        assert run_data == {}
        assert sample_data == {}


class TestBuildIndexAssignmentEdgeCases:
    """Edge cases for _build_index_assignment_from_stats."""

    def test_build_index_assignment_no_samplestats_returns_empty(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_inner, total = m._build_index_assignment_from_stats({"NumPoloniesBeforeTrimming": 1000}, "RUN01-a1b2")
        assert run_inner == {}
        assert total == 1000

    def test_build_index_assignment_sample_without_occurrences_skipped(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        # Stats with one sample that has no Occurrences
        stats = {
            "RunName": "RUN01",
            "AnalysisID": "a1b2c3d4",
            "NumPoloniesBeforeTrimming": 50000,
            "SampleStats": [
                {"SampleID": "s1", "SampleName": "S1", "NumPolonies": 100},
            ],
        }
        run_inner, total = m._build_index_assignment_from_stats(stats, "RUN01-a1b2")
        assert run_inner == {}
        assert total == 50000


class TestMergeManifestIndexSequencesEdgeCases:
    """Edge cases for _merge_manifest_index_sequences."""

    def test_merge_manifest_no_samples_returns_early(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        sample_to_index = {"RUN01-a1b2": {"AAATTT": {}}}
        m._merge_manifest_index_sequences(sample_to_index, {}, "RUN01-a1b2")
        assert sample_to_index["RUN01-a1b2"]["AAATTT"].get("Index1", "") == ""

    def test_merge_manifest_run_not_in_assignment_returns_early(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        sample_to_index = {}
        m._merge_manifest_index_sequences(
            sample_to_index,
            {"Samples": [{"SampleName": "S1", "Indexes": [{"Index1": "A", "Index2": "T"}]}]},
            "RUN01-a1b2",
        )
        assert sample_to_index == {}

    def test_merge_manifest_merged_indices_not_in_run_data_skipped(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_inner = {"AAATTT": {"SampleID": "RUN01-a1b2__S1"}}
        sample_to_index = {"RUN01-a1b2": run_inner}
        m._merge_manifest_index_sequences(
            sample_to_index,
            {"Samples": [{"SampleName": "S1", "Indexes": [{"Index1": "XXX", "Index2": "YYY"}]}]},
            "RUN01-a1b2",
        )
        assert run_inner["AAATTT"].get("Index1", "") == ""


class TestParseRunUnassignedEdgeCases:
    """Edge cases for _parse_run_unassigned_sequences."""

    def test_parse_run_unassigned_empty_data_source_returns_empty(self, fixtures_dir):
        report.reset()
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject")]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        assert m._parse_run_unassigned_sequences("") == {}

    def test_parse_run_unassigned_no_lanes_skipped(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndNoProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        unassigned = m._parse_run_unassigned_sequences("bases2fastq/run")
        assert unassigned == {}


class TestModuleNoSamplesFound:
    """Tests that ModuleNoSamplesFound is raised when no data."""

    def test_no_log_files_raises(self, tmp_path):
        report.reset()
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        report.analysis_files = [str(empty_dir)]
        report.search_files(["bases2fastq"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()


class TestProjectLevelPath:
    """Tests for project_level summary path (tabulate_project_stats, manifest in project)."""

    def test_project_level_only_produces_sections(self, fixtures_dir, tmp_path):
        """Directory with only project-level RunStats (no run-level) uses project_level path."""
        report.reset()
        project_stats = _load_fixture(
            fixtures_dir, "PairedEndDefaultProject", "Samples", "DefaultProject", "DefaultProject_RunStats.json"
        )
        manifest = _load_fixture(fixtures_dir, "PairedEndDefaultProject", "RunManifest.json")
        (tmp_path / "Samples" / "DefaultProject").mkdir(parents=True)
        (tmp_path / "Samples" / "DefaultProject" / "DefaultProject_RunStats.json").write_text(json.dumps(project_stats))
        (tmp_path / "RunManifest.json").write_text(json.dumps(manifest))
        report.analysis_files = [str(tmp_path)]
        report.search_files(["bases2fastq"])
        config.strict = True
        m = MultiqcModule()
        assert len(m.project_level_data) >= 1
        assert len(m.run_level_data) == 0
        assert len(m.sections) > 0


class TestSelectDataBySummaryPath:
    """Tests for _select_data_by_summary_path branches."""

    def test_select_data_project_level(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_data, sample_data, samples_to_projects, manifest_data, index_data, unassigned = (
            m._select_data_by_summary_path("project_level")
        )
        assert run_data is m.project_level_data
        assert sample_data is m.project_level_samples
        assert unassigned == {}

    def test_select_data_combined_level(self, fixtures_dir):
        report.reset()
        run_dir = fixtures_dir / "PairedEndDefaultProject"
        report.analysis_files = [str(run_dir)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        run_data, sample_data, samples_to_projects, manifest_data, index_data, unassigned = (
            m._select_data_by_summary_path("combined_level")
        )
        assert run_data is m.run_level_data
        assert sample_data is m.project_level_samples
        assert isinstance(unassigned, dict)


class TestParseIndexAssignmentEdgeCases:
    """Edge cases for _parse_index_assignment."""

    def test_parse_index_assignment_runstats_missing_samplestats(self, fixtures_dir, tmp_path):
        report.reset()
        run_stats = {"RunName": "RUN01", "AnalysisID": "a1b2c3d4", "NumPolonies": 100}
        (tmp_path / "RunStats.json").write_text(json.dumps(run_stats))
        (tmp_path / "RunManifest.json").write_text(json.dumps({"Settings": [{"Lane": 1}]}))
        report.analysis_files = [str(fixtures_dir / "PairedEndNoProject"), str(tmp_path)]
        report.search_files(["bases2fastq"])
        m = MultiqcModule()
        result = m._parse_index_assignment("bases2fastq/manifest")
        assert isinstance(result, dict)


def _test_data_bases2fastq_dir():
    """Path to test-data/data/modules/bases2fastq (used for skipif, no fixture)."""
    repo_root = Path(__file__).resolve().parents[4]
    return repo_root / "test-data" / "data" / "modules" / "bases2fastq"


class TestIntegration:
    """Integration test using test-data repo (skipped when absent)."""

    @pytest.mark.skipif(
        not _test_data_bases2fastq_dir().exists(),
        reason="test-data/data/modules/bases2fastq not found (clone test-data repo)",
    )
    def test_module_run_with_test_data(self, data_dir):
        """Full module run against test-data repo produces sections and general stats."""
        report.reset()
        mod_dir = data_dir / "modules" / "bases2fastq"
        report.analysis_files = [str(mod_dir)]
        report.search_files(["bases2fastq"])
        config.strict = True
        m = MultiqcModule()
        # Test-data has multiple run roots (WGS, WES, PairedEndNoProject, PairedEndDefaultProject, etc.)
        total_samples = len(m.run_level_samples) + len(m.project_level_samples)
        assert len(m.run_level_data) >= 2 or len(m.project_level_data) >= 1, (
            "expected at least 2 runs or at least 1 project from test-data"
        )
        assert total_samples >= 10, "expected at least 10 samples from test-data"
        # Module must produce output (general stats and/or sections)
        assert len(report.general_stats_data) > 0 or len(m.sections) > 0, (
            "expected general stats or report sections to be populated"
        )
