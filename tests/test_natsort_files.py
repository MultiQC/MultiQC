#!/usr/bin/env python
"""
Tests for natural sorting of files by clean sample name
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from multiqc.base_module import BaseMultiqcModule
from multiqc.types import FileDict, ModuleId


def test_files_natsorted_by_clean_sample_name():
    """Test that files are naturally sorted by clean sample name"""
    # Create mock files with sample names that should be sorted
    mock_files = [
        FileDict({"fn": "sample_10.txt", "root": "root_dir", "sp_key": "test"}),
        FileDict({"fn": "sample_2.txt", "root": "root_dir", "sp_key": "test"}),
        FileDict({"fn": "sample_1.txt", "root": "root_dir", "sp_key": "test"}),
    ]

    # Expected order after natural sorting
    expected_order = ["sample_1.txt", "sample_2.txt", "sample_10.txt"]

    # Mock report.files
    with patch("multiqc.report.files", {ModuleId("test"): mock_files}):
        # Create a test module
        test_module = BaseMultiqcModule(name="Test Module", anchor="test")

        # Get sorted files
        sorted_files = []
        for f in test_module.find_log_files("test", filecontents=False, filehandles=False):
            sorted_files.append(f["fn"])

        # Assert files are in correct order
        assert sorted_files == expected_order, (
            f"Files not sorted correctly. Expected {expected_order}, got {sorted_files}"
        )


def test_files_natsorted_with_complex_sample_names():
    """Test that files with complex naming patterns are naturally sorted by clean sample name"""
    # Create mock files with various naming patterns
    mock_files = [
        FileDict({"fn": "SRR10_1.fastq.gz", "root": "root_dir", "sp_key": "test"}),
        FileDict({"fn": "SRR2_1.fastq.gz", "root": "root_dir", "sp_key": "test"}),
        FileDict({"fn": "SRR100_1.fastq.gz", "root": "root_dir", "sp_key": "test"}),
    ]

    # Expected order after natural sorting
    expected_order = ["SRR2_1.fastq.gz", "SRR10_1.fastq.gz", "SRR100_1.fastq.gz"]

    # Mock report.files
    with patch("multiqc.report.files", {ModuleId("test"): mock_files}):
        # Create a test module
        test_module = BaseMultiqcModule(name="Test Module", anchor="test")

        # Get sorted files
        sorted_files = []
        for f in test_module.find_log_files("test", filecontents=False, filehandles=False):
            sorted_files.append(f["fn"])

        # Assert files are in correct order
        assert sorted_files == expected_order, (
            f"Files not sorted correctly. Expected {expected_order}, got {sorted_files}"
        )


if __name__ == "__main__":
    test_files_natsorted_by_clean_sample_name()
    test_files_natsorted_with_complex_sample_names()
    print("All tests passed!")
