"""Tests for the FastQE MultiQC module."""
from pathlib import Path

import pytest
from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.fastqe import MultiqcModule
SAMPLE_FASTQE_DATA = (
    "Filename\tStatistic\tQualities\n"
    "sample1.fastq.gz\tmean\t😁😁😁😁😁😄😄😄😃😃😃😉😉😊😊😊😊😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋😋\n"
)


def run_module(analysis_dir: Path, monkeypatch=None) -> MultiqcModule:
    report.reset()
    if monkeypatch is not None:
        monkeypatch.setattr(config, "preserve_module_raw_data", True)
    else:
        config.preserve_module_raw_data = True
    report.analysis_files = [str(analysis_dir)]
    report.search_files(["fastqe"])
    return MultiqcModule()


class TestFastqeParsing:
    def test_basic_parse(self, tmp_path, monkeypatch):
        (tmp_path / "sample1_fastqe.tsv").write_text(SAMPLE_FASTQE_DATA)
        m = run_module(tmp_path, monkeypatch)
        assert len(m.saved_raw_data["multiqc_fastqe"]) > 0
    def test_sample_name(self, tmp_path, monkeypatch):
        (tmp_path / "sample1_fastqe.tsv").write_text(SAMPLE_FASTQE_DATA)
        m = run_module(tmp_path, monkeypatch)
        data = m.saved_raw_data["multiqc_fastqe"]
        assert "sample1" in data

    def test_quality_string_present(self, tmp_path, monkeypatch):
        (tmp_path / "sample1_fastqe.tsv").write_text(SAMPLE_FASTQE_DATA)
        m = run_module(tmp_path, monkeypatch)
        data = m.saved_raw_data["multiqc_fastqe"]
        for s_name, stats in data.items():
            assert "mean" in stats
            assert len(stats["mean"]) > 0
    def test_section_created(self, tmp_path, monkeypatch):
        (tmp_path / "sample1_fastqe.tsv").write_text(SAMPLE_FASTQE_DATA)
        m = run_module(tmp_path, monkeypatch)
        assert len(m.sections) == 1


class TestFastqeEdgeCases:
    def test_no_matching_files(self, tmp_path):
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_empty_file(self, tmp_path):
        (tmp_path / "test_fastqe.txt").write_text("")
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_header_only(self, tmp_path):
        (tmp_path / "test_fastqe.txt").write_text("Filename\tStatistic\tQualities\n")
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()


class TestFastqeMultipleSamples:
    def test_multi_sample_file(self, tmp_path, monkeypatch):
        test_data = (
            "Filename\tStatistic\tQualities\n"
            "sample1.fq\tmean\t😘 😆 😌 😌 😋 😜 😜 😜\n"
            "sample2.fq\tmean\t😋 😜 😜 😆 😘 😌 😌 😜\n"
        )
        (tmp_path / "multi_fastqe.txt").write_text(test_data)
        m = run_module(tmp_path, monkeypatch)
        data = m.saved_raw_data["multiqc_fastqe"]
        assert len(data) == 2
        assert "sample1" in data
        assert "sample2" in data
