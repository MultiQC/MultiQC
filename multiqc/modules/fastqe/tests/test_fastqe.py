# Test data generated with:
# !pip install fastqe
# !wget https://raw.githubusercontent.com/fastqe/fastqe/refs/heads/master/tests/data/test.fastq
# !fastqe test.fastq --output test_fastqe.txt

"""
Filename        Statistic       Qualities
test.fastq      mean    😘 😆 😌 😌 😋 😜 😜 😜 😋 😜 😜 😉 😁 😁 😁 😁 😁 😉 😉 😁 😁 😁 😄 😁 😁 😄 😄 😁 😄 😁 😁 😁 😄 😄 😄 😄 😁 😄 😉 😁 😄 😁 😁 😄 😄 😁 😁 😄 😉 😁 😁 😝 😝 😁 😁 😌 😁 😉 😁 😁 😁 😁 😁 😝 😁 😁 😁 😁 😝 😁 😁 😉 😁 😁 😉 😁 😉 😃 😉 😁 😁 😁 😄 😁 😄 😄 😁 😄 😁 😝 😁 😄 😄 😁 😉 😁 😉 😁 😁 😁 😁 😄 😄 😁 😁 😌 😉 😜 😁 😁 😁 😁 😁 😉 😁 😜 😁 😝 😁 😁 😁 😁 😁 😁 😁 😁 😁 😉 😄 😜 😁 😄 😁 😄 😁 😄 😁 😁 😜 😜 😜 😜 😉 😜 😉 😁 😉 😋 😜 😛 😉 😜 😜 😜 😉 😉 😉 😁 😁 😁 😁 😁 😋 😛 😜 😜 😁 😛 😄 😋 😛 😝 😉 😉 😄 😌 😜 😜 😛 😆 😛 😉 😜 😜 😜 😜 😋 😝 😋 😜 😉 😉 😄 😉 😝 😝 😛 😋 😜 😜 😜 😄 😆 😜 😝 😉 😜 😋 😉 😜 😄 😌 😌 😜 😘 😝 😆 😄 😜 😜 😜 😉 😛 😄 😌 😙 😊 😄 😝 😆 😜 😉 😜 😄 😄 😌 😌 😜 😉 😜 😉 😆 😆 😛 😙 😃 😙 😘 😝 😙 😡
"""
from pathlib import Path

import pytest

from multiqc import report
from multiqc.modules.fastqe import MultiqcModule


def test_parse_fastqe(tmp_path):
    """Test parsing FastQE output"""

    # Create test input file
    test_data = """\
Filename\tStatistic\tQualities
test.fastq\tmean\t😘 😆 😌 😌 😋 😜 😜 😜 😋 😜 😜 😉
"""
    f = tmp_path / "test_fastqe.txt"
    f.write_text(test_data)

    # Run module
    report.analysis_files = [f]
    report.search_files(["fastqe"])
    m = MultiqcModule()

    # Check parsed data
    assert len(m.saved_raw_data) == 1
    data = m.saved_raw_data["multiqc_fastqe"]
    assert "test.fastq" in data
    assert "mean" in data["test.fastq"]
    assert data["test.fastq"]["mean"].startswith("😘 😆 😌")

    # Check plot data
    assert len(m.sections) == 1
    plot = m.sections[0].plot
    assert plot["id"] == "fastqe_quality_plot"
    assert plot["plot_type"] == "html"
    assert len(plot["data"]) == 1
    assert plot["data"][0]["name"] == "test.fastq"
    assert "😘 😆 😌" in plot["data"][0]["qualities"]


def test_parse_empty_file(tmp_path):
    """Test handling empty input file"""
    f = tmp_path / "empty.txt"
    f.write_text("")

    report.analysis_files = [f]
    report.search_files(["fastqe"])

    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()


def test_parse_malformed_file(tmp_path):
    """Test handling malformed input file"""
    f = tmp_path / "malformed.txt"
    f.write_text("Invalid data\nMore invalid data")

    report.analysis_files = [f]
    report.search_files(["fastqe"])

    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()


def test_multiple_samples(tmp_path):
    """Test parsing multiple samples"""
    test_data = """\
Filename\tStatistic\tQualities
sample1.fq\tmean\t😘 😆 😌
sample2.fq\tmean\t😋 😜 😜
"""
    f = tmp_path / "multi_samples.txt"
    f.write_text(test_data)

    report.analysis_files = [f]
    report.search_files(["fastqe"])
    m = MultiqcModule()

    data = m.saved_raw_data["multiqc_fastqe"]
    assert len(data) == 2
    assert "sample1.fq" in data
    assert "sample2.fq" in data

    plot = m.sections[0].plot
    assert len(plot["data"]) == 2
