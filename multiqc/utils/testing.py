import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import json
import copy
import pytest

from multiqc import config, report


def data_dir():
    test_data_dir = Path(os.environ.get("MULTIQC_TEST_DATA_DIR", config.REPO_DIR / "test-data"))
    if not test_data_dir.exists():
        raise FileNotFoundError(
            f"The test data directory expected to be found at {test_data_dir}. Please, "
            f"clone the repository with the test data by changing to the MultiQC repo root, "
            f"and running: `git clone https://github.com/MultiQC/test-data`"
        )

    return test_data_dir / "data"


class ModuleSnapshot:
    """
    A utility class for creating comprehensive snapshots of MultiQC module outputs.

    This class captures all the important data structures that modules produce,
    ensuring consistency and reproducibility in testing.
    """

    def __init__(self, module_instance, normalize_paths: bool = True, exclude_keys: Optional[List[str]] = None):
        """
        Initialize the snapshot utility.

        Args:
            module_instance: The MultiQC module instance to snapshot
            normalize_paths: Whether to normalize file paths for cross-platform compatibility
            exclude_keys: List of keys to exclude from snapshots (e.g., timestamps, temp paths)
        """
        self.module = module_instance
        self.normalize_paths = normalize_paths
        self.exclude_keys = exclude_keys or []

    def _normalize_data(self, data: Any) -> Any:
        """
        Normalize data for consistent snapshots across different environments.

        This includes:
        - Converting Path objects to strings
        - Normalizing file paths
        - Removing excluded keys
        - Sorting dictionaries for consistent ordering
        """
        if isinstance(data, dict):
            normalized = {}
            for key, value in data.items():
                if key in self.exclude_keys:
                    continue
                normalized[key] = self._normalize_data(value)
            # Sort dictionary keys for consistent ordering
            return dict(sorted(normalized.items()))
        elif isinstance(data, list):
            return [self._normalize_data(item) for item in data]
        elif isinstance(data, Path):
            path_str = str(data)
            if self.normalize_paths:
                # Normalize path separators for cross-platform compatibility
                return path_str.replace("\\", "/")
            return path_str
        elif isinstance(data, float):
            # Round floats to avoid precision issues
            return round(data, 6)
        else:
            return data

    def get_saved_raw_data(self) -> Dict[str, Any]:
        """Get the module's saved raw data, normalized for snapshotting."""
        if self.module.saved_raw_data is None:
            return {}
        return self._normalize_data(self.module.saved_raw_data)

    def get_general_stats_data(self) -> Dict[str, Any]:
        """Get general stats data from the report, normalized for snapshotting."""
        # Find this module's general stats data
        module_general_stats = {}
        for section_key, data in report.general_stats_data.items():
            if str(section_key).startswith(str(self.module.anchor)):
                module_general_stats[str(section_key)] = self._normalize_data(data)
        return module_general_stats

    def get_general_stats_headers(self) -> Dict[str, Any]:
        """Get general stats headers from the report, normalized for snapshotting."""
        module_headers = {}
        for section_key, headers in report.general_stats_headers.items():
            if str(section_key).startswith(str(self.module.anchor)):
                module_headers[str(section_key)] = self._normalize_data(headers)
        return module_headers

    def get_sections_data(self) -> List[Dict[str, Any]]:
        """Get all section data from the module, normalized for snapshotting."""
        sections_data = []
        for section in self.module.sections:
            section_data = {
                "name": section.name,
                "anchor": str(section.anchor),
                "description": section.description,
                "comment": section.comment,
                "helptext": section.helptext,
                "plot_type": section.plot.plot_type if section.plot else None,
                "plot_id": str(section.plot.id) if section.plot else None,
            }

            # Include plot data if available
            if section.plot and hasattr(section.plot, "data"):
                section_data["plot_data"] = self._normalize_data(section.plot.data)

            # Include plot configuration
            if section.plot and hasattr(section.plot, "pconfig"):
                section_data["plot_config"] = self._normalize_data(section.plot.pconfig)

            sections_data.append(section_data)

        return sections_data

    def get_software_versions(self) -> Dict[str, Any]:
        """Get software versions detected by the module."""
        return self._normalize_data(dict(self.module.versions))

    def get_complete_snapshot(self) -> Dict[str, Any]:
        """
        Get a complete snapshot of all module data.

        Returns a dictionary containing:
        - saved_raw_data: The module's parsed raw data
        - general_stats_data: Data added to the general statistics table
        - general_stats_headers: Headers for the general statistics table
        - sections: All sections added by the module
        - software_versions: Software versions detected
        - module_info: Basic module information
        """
        return {
            "module_info": {
                "name": self.module.name,
                "anchor": str(self.module.anchor),
                "info": self.module.info,
                "href": self.module.href,
                "doi": self.module.doi,
            },
            "saved_raw_data": self.get_saved_raw_data(),
            "general_stats_data": self.get_general_stats_data(),
            "general_stats_headers": self.get_general_stats_headers(),
            "sections": self.get_sections_data(),
            "software_versions": self.get_software_versions(),
        }


def run_module_test(
    module_class,
    data_files: Union[str, Path, List[Union[str, Path]]],
    module_kwargs: Optional[Dict[str, Any]] = None,
    config_updates: Optional[Dict[str, Any]] = None,
) -> ModuleSnapshot:
    """
    Run a module test and return a snapshot utility.

    Args:
        module_class: The MultiQC module class to test
        data_files: Path(s) to test data files
        module_kwargs: Additional kwargs to pass to the module constructor
        config_updates: Configuration updates to apply before running the module

    Returns:
        ModuleSnapshot instance for the tested module
    """
    # Reset report state
    report.reset()

    # Apply config updates if provided
    if config_updates:
        for key, value in config_updates.items():
            setattr(config, key, value)

    # Ensure data_files is a list
    if not isinstance(data_files, list):
        data_files = [data_files]

    # Set up analysis files
    report.analysis_files = [Path(f) for f in data_files]

    # Search for files (assuming module name matches class name)
    module_name = module_class.__module__.split(".")[-2]  # Extract module name from path
    report.search_files([module_name])

    # Enable raw data preservation for testing
    config.preserve_module_raw_data = True

    # Initialize and run the module
    module_kwargs = module_kwargs or {}
    module_instance = module_class(**module_kwargs)

    return ModuleSnapshot(module_instance)


def assert_module_data_integrity(snapshot: ModuleSnapshot):
    """
    Assert basic data integrity for a module snapshot.

    This performs common checks that should pass for any well-formed module:
    - Module has a name and anchor
    - If raw data exists, it's not empty
    - If sections exist, they have required fields
    - General stats data is properly formatted
    """
    # Check basic module info
    assert snapshot.module.name, "Module must have a name"
    assert snapshot.module.anchor, "Module must have an anchor"

    # Check raw data integrity
    raw_data = snapshot.get_saved_raw_data()
    if raw_data:
        assert isinstance(raw_data, dict), "Raw data must be a dictionary"
        for file_key, file_data in raw_data.items():
            assert isinstance(file_data, dict), f"Raw data for {file_key} must be a dictionary"
            assert len(file_data) > 0, f"Raw data for {file_key} must not be empty"

    # Check sections integrity
    sections = snapshot.get_sections_data()
    for section in sections:
        assert section["name"], "Section must have a name"
        assert section["anchor"], "Section must have an anchor"

    # Check general stats integrity
    general_stats = snapshot.get_general_stats_data()
    for section_key, data in general_stats.items():
        assert isinstance(data, dict), f"General stats data for {section_key} must be a dictionary"
        for sample_name, sample_data in data.items():
            assert isinstance(sample_data, dict), f"Sample data for {sample_name} must be a dictionary"


class BaseModuleTest:
    """
    Base class for module snapshot tests.
    
    This class provides common functionality for testing MultiQC modules
    with snapshot testing capabilities.
    """
    
    # Override these in subclasses
    MODULE_CLASS = None
    MODULE_NAME = None
    DATA_DIR_NAME = None  # If different from MODULE_NAME
    
    def get_module_data_dir(self, data_dir):
        """Get the data directory for this module."""
        dir_name = self.DATA_DIR_NAME or self.MODULE_NAME
        return data_dir / "modules" / dir_name
    
    def create_module_snapshot(self, data_dir, snapshot_config):
        """
        Create a snapshot of the module with all available test data.
        
        This method automatically finds all test files in the module's data directory
        and runs the module on them, returning a ModuleSnapshot instance.
        """
        if self.MODULE_CLASS is None:
            pytest.skip("MODULE_CLASS not defined")
        
        module_data_dir = self.get_module_data_dir(data_dir)
        
        # Apply snapshot configuration
        from multiqc import config
        for key, value in snapshot_config.items():
            setattr(config, key, value)
        
        # Run the module test
        return run_module_test(
            module_class=self.MODULE_CLASS,
            data_files=list(module_data_dir.rglob("*")),
            config_updates=snapshot_config,
        )
    
    def assert_module_data_integrity(self, module_snapshot):
        """Test basic data integrity for the module."""
        assert_module_data_integrity(module_snapshot)
