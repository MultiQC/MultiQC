"""Test that the sample renaming functionality works as expected."""

import tempfile
import os
import subprocess
import json
from pathlib import Path

import pytest

from multiqc import config
from multiqc.core.update_config import ClConfig, update_config


@pytest.fixture
def sample_names_file():
    """Create a temporary sample names file for testing."""
    content = """MultiQC Names\tOriginal Sample Name\tid_Original Sample Name
Sample_A\tOriginal_A\tID_A
Sample_B\tOriginal_B\tID_B
Sample_C\tOriginal_C\tID_C"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        f.flush()
        yield f.name
    os.unlink(f.name)


def test_load_sample_names_file(sample_names_file):
    """Test that sample names file is loaded correctly."""
    # Reset config
    config.reset()

    # Load the sample names file
    config.load_sample_names(Path(sample_names_file))

    # Check that the config was updated correctly
    assert len(config.sample_names_rename_buttons) == 3
    assert config.sample_names_rename_buttons == ["MultiQC Names", "Original Sample Name", "id_Original Sample Name"]

    assert len(config.sample_names_rename) == 3
    assert config.sample_names_rename[0] == ["Sample_A", "Original_A", "ID_A"]
    assert config.sample_names_rename[1] == ["Sample_B", "Original_B", "ID_B"]
    assert config.sample_names_rename[2] == ["Sample_C", "Original_C", "ID_C"]


def test_sample_names_config_integration():
    """Test that sample names config is properly integrated via update_config."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("Button1\tButton2\nSample1\tAltName1\nSample2\tAltName2")
        f.flush()

        try:
            # Reset config
            config.reset()

            # Use update_config to load sample names
            update_config(cfg=ClConfig(sample_names=f.name))

            # Verify the config was updated
            assert len(config.sample_names_rename_buttons) == 2
            assert config.sample_names_rename_buttons == ["Button1", "Button2"]
            assert len(config.sample_names_rename) == 2
            assert config.sample_names_rename[0] == ["Sample1", "AltName1"]
            assert config.sample_names_rename[1] == ["Sample2", "AltName2"]
        finally:
            os.unlink(f.name)


def test_sample_names_cmdline_integration(tmp_path):
    """Test that sample renaming works via command line interface with custom content."""

    # Create custom content with known sample names
    custom_content = "#id: test_custom\n#plot_type: table\n#section_name: Test Custom Data\nSample_Name\tValue1\tValue2\nSample1\t100\t200\nSample2\t150\t250\nSample3\t120\t220"

    custom_file = tmp_path / "test_data_mqc.txt"
    with open(custom_file, "w") as f:
        f.write(custom_content)

    # Create sample names file that maps the actual sample names
    sample_names_content = "MultiQC Names\tOriginal Sample Name\tid_Original Sample Name\nSample1\tOriginal_Sample_1\tID_Sample_1\nSample2\tOriginal_Sample_2\tID_Sample_2\nSample3\tOriginal_Sample_3\tID_Sample_3"

    sample_names_file = tmp_path / "sample_names.txt"
    with open(sample_names_file, "w") as f:
        f.write(sample_names_content)

    # Run MultiQC with sample names file
    cmd = ["multiqc", str(tmp_path), "--sample-names", str(sample_names_file), "--filename", "test_report", "--force"]

    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)

    # Check that MultiQC ran successfully
    assert result.returncode == 0, f"MultiQC failed with: {result.stderr}"

    # Check that report was generated
    assert (tmp_path / "test_report.html").exists()
    assert (tmp_path / "test_report_data").exists()

    # Check that the data file contains the sample rename configuration
    with open(tmp_path / "test_report_data" / "multiqc_data.json", "r") as f:
        data = json.load(f)

    # Verify that sample_names_rename config is present in the report data
    assert "config_sample_names_rename" in data
    assert "config_sample_names_rename_buttons" in data

    # Verify the content matches what we expect
    assert data["config_sample_names_rename_buttons"] == [
        "MultiQC Names",
        "Original Sample Name",
        "id_Original Sample Name",
    ]
    assert len(data["config_sample_names_rename"]) == 3

    # Verify actual mapping content
    assert data["config_sample_names_rename"][0] == ["Sample1", "Original_Sample_1", "ID_Sample_1"]
    assert data["config_sample_names_rename"][1] == ["Sample2", "Original_Sample_2", "ID_Sample_2"]
    assert data["config_sample_names_rename"][2] == ["Sample3", "Original_Sample_3", "ID_Sample_3"]

    # Verify that HTML report contains sample renaming buttons
    with open(tmp_path / "test_report.html", "r") as f:
        html_content = f.read()

    # Check that sample rename buttons are present
    assert "Change sample names:" in html_content
    assert "MultiQC Names" in html_content
    assert "Original Sample Name" in html_content
    assert "id_Original Sample Name" in html_content
    assert "mqc_sname_switches" in html_content

    # Check that the initial state shows the correct sample names in the table content
    # When "MultiQC Names" is active (index 0), it should show Sample1, Sample2, Sample3
    import re

    # Extract the table content (the actual displayed data)
    table_matches = re.findall(r'<span class="th-sample-name"[^>]*>([^<]+)</span>', html_content)
    assert len(table_matches) == 3, f"Expected 3 sample names in table, got {len(table_matches)}: {table_matches}"

    # Verify that the table shows the original sample names (Sample1, Sample2, Sample3)
    assert "Sample1" in table_matches
    assert "Sample2" in table_matches
    assert "Sample3" in table_matches

    # Verify that the renamed versions are NOT displayed in the table content
    assert "Original_Sample_1" not in table_matches
    assert "ID_Sample_1" not in table_matches

    # The config data should still contain the mapping (this is expected for the buttons to work)
    assert '"sample_names_rename"' in html_content
    assert '"Original_Sample_1"' in html_content  # This should be in the config, not the table


def test_sample_names_malformed_file():
    """Test handling of malformed sample names files."""

    # Test file with inconsistent columns
    malformed_content = """Header1\tHeader2\tHeader3
Sample1\tAlt1\tID1
Sample2\tAlt2
Sample3\tAlt3\tID3\tExtra"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(malformed_content)
        f.flush()

        try:
            config.reset()
            config.load_sample_names(Path(f.name))

            # Should still load header and all lines (but warn about inconsistent columns)
            assert len(config.sample_names_rename_buttons) == 3
            # All lines should be loaded but some may be malformed
            assert len(config.sample_names_rename) == 3  # All three lines loaded, despite warnings

        finally:
            os.unlink(f.name)


def test_sample_names_empty_file():
    """Test handling of empty sample names file."""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("")
        f.flush()

        try:
            config.reset()
            config.load_sample_names(Path(f.name))

            # Should have empty rename configs
            assert len(config.sample_names_rename_buttons) == 0
            assert len(config.sample_names_rename) == 0

        finally:
            os.unlink(f.name)


def test_sample_names_single_column_file():
    """Test that files with single columns are handled correctly."""

    content = """OnlyHeader
Sample1
Sample2"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        f.flush()

        try:
            config.reset()
            config.load_sample_names(Path(f.name))

            # Single column files shouldn't be processed as rename configs
            assert len(config.sample_names_rename_buttons) == 0
            assert len(config.sample_names_rename) == 0

        finally:
            os.unlink(f.name)
