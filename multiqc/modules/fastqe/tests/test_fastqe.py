"""Tests for the FastQE MultiQC module."""

from pathlib import Path

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.fastqe import MultiqcModule


@pytest.fixture()
def data_dir():
    return Path(__file__).parent / "data"


def run_module(analysis_dir: Path) -> MultiqcModule:
    report.reset()
    config.preserve_module_raw_data = True
    report.analysis_files = [str(analysis_dir)]
    report.search_files(["fastqe"])
    return MultiqcModule()


class TestFastqeParsing:
    def test_basic_parse(self, data_dir):
        m = run_module(data_dir)
        assert len(m.saved_raw_data["multiqc_fastqe"]) > 0

    def test_sample_name(self, data_dir):
        m = run_module(data_dir)
        data = m.saved_raw_data["multiqc_fastqe"]
        found_samples = list(data.keys())
        assert len(found_samples) >= 1

    def test_quality_string_present(self, data_dir):
        m = run_module(data_dir)
        data = m.saved_raw_data["multiqc_fastqe"]
        for s_name, stats in data.items():
            assert "mean" in stats
            assert len(stats["mean"]) > 0

    def test_sections_created(self, data_dir):
        m = run_module(data_dir)
        assert len(m.sections) >= 1

    def test_general_stats(self, data_dir):
        run_module(data_dir)
        assert len(report.general_stats_data) > 0

    def test_numeric_stats(self, data_dir):
        m = run_module(data_dir)
        assert len(m.sections) == 2


class TestFastqeEdgeCases:
    def test_no_matching_files(self, tmp_path):
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_empty_file(self, tmp_path):
        f = tmp_path / "test_fastqe.txt"
        f.write_text("")
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    def test_header_only(self, tmp_path):
        f = tmp_path / "test_fastqe.txt"
        f.write_text("Filename\tStatistic\tQualities\n")
        report.reset()
        report.analysis_files = [str(tmp_path)]
        report.search_files(["fastqe"])
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()


class TestFastqeMultipleSamples:
    def test_multi_sample_file(self, tmp_path):
        test_data = (
            "Filename\tStatistic\tQualities\n"
            "sample1.fq\tmean\t😘 😆 😌 😌 😋 😜 😜 😜\n"
            "sample2.fq\tmean\t😋 😜 😜 😆 😘 😌 😌 😜\n"
        )
        f = tmp_path / "multi_fastqe.txt"
        f.write_text(test_data)
        m = run_module(tmp_path)

        data = m.saved_raw_data["multiqc_fastqe"]
        assert len(data) == 2
        # clean_s_name strips the .fq extension
        assert "sample1" in data
        assert "sample2" in data
