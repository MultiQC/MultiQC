"""
Template for MultiQC module snapshot tests.

This file serves as a template for implementing comprehensive snapshot tests
for MultiQC modules. Copy this file and modify it for your specific module.

To use this template:
1. Copy this file to your module's test directory
2. Replace 'YourModule' with your actual module name
3. Update MODULE_CLASS and MODULE_NAME
4. Customize the test methods as needed
5. Run the tests to generate initial snapshots
"""

import pytest
from pathlib import Path

from multiqc import config
from multiqc.utils import testing
from tests.conftest import BaseModuleTest

# Import your module class here
# from multiqc.modules.yourmodule import MultiqcModule


class TestYourModuleSnapshot(BaseModuleTest):
    """
    Comprehensive snapshot tests for YourModule.

    This class provides complete snapshot testing for a MultiQC module,
    ensuring that all outputs are documented and reproducible.
    """

    # Update these for your module
    MODULE_CLASS = None  # Replace with your MultiqcModule class
    MODULE_NAME = "yourmodule"  # Replace with your module name
    DATA_DIR_NAME = None  # Optional: if data dir name differs from module name

    def test_module_basic_functionality(self, module_snapshot):
        """Test basic module functionality and data integrity."""
        # Basic integrity checks
        testing.assert_module_data_integrity(module_snapshot)

        # Module-specific checks
        raw_data = module_snapshot.get_saved_raw_data()

        # Example: Check that we have the expected data files
        # assert "multiqc_yourmodule" in raw_data

        # Example: Check that we have data for expected samples
        # expected_samples = {"sample1", "sample2"}
        # actual_samples = set(raw_data["multiqc_yourmodule"].keys())
        # assert expected_samples.issubset(actual_samples)

    def test_module_parsing_consistency(self, module_snapshot):
        """Test that parsing is consistent across different input files."""
        raw_data = module_snapshot.get_saved_raw_data()

        # Example: Check that all samples have consistent data structure
        # if "multiqc_yourmodule" in raw_data:
        #     sample_data = raw_data["multiqc_yourmodule"]
        #     if len(sample_data) > 1:
        #         sample_keys = [set(data.keys()) for data in sample_data.values()]
        #         first_keys = sample_keys[0]
        #         for keys in sample_keys[1:]:
        #             assert keys == first_keys, "Inconsistent parsing across samples"

    def test_module_data_snapshot(self, module_snapshot, snapshot):
        """Snapshot test for module raw data."""
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data == snapshot

    def test_module_general_stats_snapshot(self, module_snapshot, snapshot):
        """Snapshot test for general statistics."""
        general_stats = {
            "data": module_snapshot.get_general_stats_data(),
            "headers": module_snapshot.get_general_stats_headers(),
        }
        assert general_stats == snapshot

    def test_module_sections_snapshot(self, module_snapshot, snapshot):
        """Snapshot test for module sections."""
        sections = module_snapshot.get_sections_data()
        assert sections == snapshot

    def test_module_complete_snapshot(self, module_snapshot, snapshot):
        """Complete snapshot test for the entire module."""
        complete_snapshot = module_snapshot.get_complete_snapshot()
        assert complete_snapshot == snapshot


class TestYourModuleSpecificFiles:
    """
    Test specific files or scenarios for YourModule.

    Use this class for testing specific input files, edge cases,
    or particular configurations.
    """

    @pytest.fixture
    def module_data_dir(self, data_dir):
        """Get the data directory for this module."""
        return data_dir / "modules" / "yourmodule"  # Update with your module name

    @pytest.mark.parametrize(
        "filename",
        [
            # Add specific test files here
            # "sample1.log",
            # "sample2.txt",
        ],
    )
    def test_individual_file(self, module_data_dir, filename, snapshot, snapshot_config):
        """Test parsing of individual files."""
        file_path = module_data_dir / filename

        if not file_path.exists():
            pytest.skip(f"Test file {filename} not found")

        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Run the module test on just this file
        module_snapshot = testing.run_module_test(
            module_class=None,  # Replace with your MODULE_CLASS
            data_files=[file_path],
            config_updates=snapshot_config,
        )

        # Create a snapshot of this file's data
        file_snapshot = {
            "filename": filename,
            "parsed_data": module_snapshot.get_saved_raw_data(),
            "general_stats": module_snapshot.get_general_stats_data(),
            "sections": module_snapshot.get_sections_data(),
        }

        assert file_snapshot == snapshot

    def test_empty_input(self, snapshot_config):
        """Test module behavior with no input files."""
        from multiqc.base_module import ModuleNoSamplesFound

        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Test that module raises ModuleNoSamplesFound with no input
        with pytest.raises(ModuleNoSamplesFound):
            testing.run_module_test(
                module_class=None,  # Replace with your MODULE_CLASS
                data_files=[],
                config_updates=snapshot_config,
            )

    def test_module_with_custom_config(self, module_data_dir, snapshot, snapshot_config):
        """Test module with custom configuration."""
        # Example: Test with custom configuration
        custom_config = {
            **snapshot_config,
            # Add your custom config here
            # "yourmodule_custom_option": True,
        }

        # Apply configuration
        for key, value in custom_config.items():
            setattr(config, key, value)

        # Run the module test
        module_snapshot = testing.run_module_test(
            module_class=None,  # Replace with your MODULE_CLASS
            data_files=list(module_data_dir.rglob("*")),
            config_updates=custom_config,
        )

        # Snapshot the results
        custom_snapshot = module_snapshot.get_complete_snapshot()
        assert custom_snapshot == snapshot


# Example of testing specific parsing functions
class TestYourModuleParsingFunctions:
    """
    Test specific parsing functions in isolation.

    This is useful for testing parsing logic without running the full module.
    """

    def test_parse_function_example(self, snapshot):
        """Test a specific parsing function."""
        # Example test data
        test_input = """
        # Example log file content
        Sample: test_sample
        Value1: 100
        Value2: 200
        """

        # Import and test your parsing function
        # from multiqc.modules.yourmodule.yourmodule import parse_log_file
        # result = parse_log_file(test_input)
        # assert result == snapshot

        # Placeholder for template
        result = {"placeholder": "Replace with actual parsing test"}
        assert result == snapshot


# Example of testing table format consistency
class TestYourModuleTableFormat:
    """
    Test table format consistency for YourModule.

    This ensures that table outputs follow consistent formatting rules.
    """

    def test_table_headers_consistency(self, module_data_dir, snapshot_config):
        """Test that table headers are consistent and well-formatted."""
        # Apply configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Run the module test
        module_snapshot = testing.run_module_test(
            module_class=None,  # Replace with your MODULE_CLASS
            data_files=list(module_data_dir.rglob("*")),
            config_updates=snapshot_config,
        )

        # Check general stats headers
        headers = module_snapshot.get_general_stats_headers()

        for section_key, section_headers in headers.items():
            for header_key, header_config in section_headers.items():
                # Check that required fields are present
                assert "title" in header_config, f"Header {header_key} missing title"
                assert "description" in header_config, f"Header {header_key} missing description"

                # Check that title is properly formatted
                title = header_config["title"]
                assert isinstance(title, str), f"Header {header_key} title must be string"
                assert len(title) > 0, f"Header {header_key} title cannot be empty"

                # Add more header validation as needed

    def test_data_format_consistency(self, module_data_dir, snapshot_config):
        """Test that data formats are consistent across samples."""
        # Apply configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Run the module test
        module_snapshot = testing.run_module_test(
            module_class=None,  # Replace with your MODULE_CLASS
            data_files=list(module_data_dir.rglob("*")),
            config_updates=snapshot_config,
        )

        # Check data consistency
        raw_data = module_snapshot.get_saved_raw_data()
        general_stats = module_snapshot.get_general_stats_data()

        # Example: Check that numeric values are properly formatted
        for file_key, file_data in raw_data.items():
            for sample_name, sample_data in file_data.items():
                for metric_name, metric_value in sample_data.items():
                    if isinstance(metric_value, (int, float)):
                        # Check that numeric values are reasonable
                        assert not (isinstance(metric_value, float) and (metric_value != metric_value)), (
                            f"Invalid numeric value for {sample_name}.{metric_name}: {metric_value}"
                        )
