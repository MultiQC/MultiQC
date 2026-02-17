import pytest

from multiqc.modules.samtools.coverage import parse_single_report


@pytest.fixture
def valid_coverage_report():
    """A valid samtools coverage report."""
    return {
        "f": (
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t100\t50\t80\t80.0\t5.5\t30.0\t40.0\n"
            "chr2\t1\t200\t100\t150\t75.0\t4.2\t28.0\t35.0\n"
        ),
        "s_name": "test_sample",
    }


@pytest.fixture
def report_with_extra_fields():
    """A samtools coverage report with extra fields (e.g., user-added 'all' line)."""
    return {
        "f": (
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t100\t50\t80\t80.0\t5.5\t30.0\t40.0\n"
            "all\t1\t300\t150\t230\t76.67\t4.85\t29.0\t37.5\textra_field\n"
        ),
        "s_name": "test_sample",
    }


@pytest.fixture
def report_with_fewer_fields():
    """A samtools coverage report with missing fields."""
    return {
        "f": (
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t100\t50\t80\t80.0\t5.5\t30.0\t40.0\n"
            "chr2\t1\t200\t100\t150\t75.0\t4.2\n"
        ),
        "s_name": "test_sample",
    }


def test_valid_coverage_report(valid_coverage_report):
    """Test that valid coverage reports are parsed correctly."""
    result = parse_single_report(valid_coverage_report)

    assert len(result) == 2
    assert "chr1" in result
    assert "chr2" in result
    assert result["chr1"]["startpos"] == 1
    assert result["chr1"]["endpos"] == 100
    assert result["chr1"]["numreads"] == 50
    assert result["chr1"]["coverage"] == 80.0
    assert result["chr2"]["numreads"] == 100


def test_report_with_extra_fields(report_with_extra_fields):
    """Test that lines with extra fields are skipped without crashing.

    This is a regression test for issue #3343 where extra fields would cause
    ValueError: too many values to unpack.
    """
    result = parse_single_report(report_with_extra_fields)

    # Only the valid line should be parsed
    assert len(result) == 1
    assert "chr1" in result
    # The 'all' line with extra fields should have been skipped
    assert "all" not in result


def test_report_with_fewer_fields(report_with_fewer_fields):
    """Test that lines with fewer fields are skipped without crashing."""
    result = parse_single_report(report_with_fewer_fields)

    # Only the valid line should be parsed
    assert len(result) == 1
    assert "chr1" in result
    # The chr2 line with fewer fields should have been skipped
    assert "chr2" not in result
