"""
Tests for the fasterp module.

Since fasterp has the same interface as fastp, these tests verify that
the fasterp module correctly processes fasterp JSON output files.
"""

import json
from pathlib import Path

import pytest

from multiqc import report
from multiqc.modules.fasterp import MultiqcModule


def test_fasterp_module(data_dir, tmp_path):
    """Test that the fasterp module can parse fasterp JSON output."""
    # Setup test data directory
    fasterp_dir = data_dir / "modules" / "fasterp"
    assert fasterp_dir.exists(), f"Test data directory not found: {fasterp_dir}"

    # Reset report state
    report.reset()
    report.analysis_files = [fasterp_dir]
    report.search_files(["fasterp"])

    # Initialize the fasterp module
    module = MultiqcModule()

    # Verify module found samples
    assert len(module.fasterp_data) > 0, "No samples found in fasterp data"

    # Verify general stats were added
    assert len(report.general_stats_data) > 0, "No general stats data found"

    # Verify sections were added
    assert len(module.sections) > 0, "No sections added to module"

    # Expected sections
    expected_section_names = [
        "Filtered Reads",
        "Duplication Rates",
        "Insert Sizes",
        "Sequence Quality",
        "GC Content",
        "N content",
        "Overrepresented Sequences",
    ]

    section_names = [s.name for s in module.sections]
    for expected in expected_section_names:
        assert expected in section_names, f"Expected section '{expected}' not found"


def test_fasterp_json_parsing(data_dir):
    """Test that fasterp correctly parses the JSON format."""
    fasterp_dir = data_dir / "modules" / "fasterp"
    json_file = fasterp_dir / "sample1.json"

    assert json_file.exists(), f"Test JSON file not found: {json_file}"

    # Verify the JSON contains the fasterp command
    with open(json_file) as f:
        data = json.load(f)

    assert "command" in data, "JSON missing 'command' field"
    assert "fasterp" in data["command"], "Command should contain 'fasterp'"
    assert "before_filtering" in data.get("summary", {}), "JSON missing 'before_filtering' data"
    assert "after_filtering" in data.get("summary", {}), "JSON missing 'after_filtering' data"


def test_fasterp_data_structure(data_dir):
    """Test that fasterp module extracts expected data fields."""
    report.reset()
    report.analysis_files = [data_dir / "modules" / "fasterp"]
    report.search_files(["fasterp"])

    module = MultiqcModule()

    # Get first sample name
    sample_name = list(module.fasterp_data.keys())[0]
    sample_data = module.fasterp_data[sample_name]

    # Verify expected data fields are present
    expected_fields = [
        "before_filtering_total_reads",
        "after_filtering_total_reads",
        "filtering_result_passed_filter_reads",
        "pct_surviving",
    ]

    for field in expected_fields:
        assert field in sample_data, f"Expected field '{field}' not found in sample data"
        assert isinstance(sample_data[field], (int, float)), f"Field '{field}' should be numeric"


def test_fasterp_vs_fastp_compatibility(data_dir):
    """Test that fasterp and fastp modules handle the same JSON format."""
    from multiqc.modules.fastp import MultiqcModule as FastpModule

    # Test with fasterp data
    report.reset()
    report.analysis_files = [data_dir / "modules" / "fasterp"]
    report.search_files(["fasterp"])
    fasterp_module = MultiqcModule()

    # Get fasterp sample data structure
    fasterp_sample = list(fasterp_module.fasterp_data.keys())[0]
    fasterp_keys = set(fasterp_module.fasterp_data[fasterp_sample].keys())

    # Test with fastp data
    report.reset()
    report.analysis_files = [data_dir / "modules" / "fastp"]
    report.search_files(["fastp"])
    fastp_module = FastpModule()

    # Get fastp sample data structure
    fastp_sample = list(fastp_module.fastp_data.keys())[0]
    fastp_keys = set(fastp_module.fastp_data[fastp_sample].keys())

    # Verify both modules extract similar data fields
    # They should have significant overlap in the data structure
    common_keys = fasterp_keys & fastp_keys
    assert len(common_keys) > 5, "fasterp and fastp should extract similar data fields"
