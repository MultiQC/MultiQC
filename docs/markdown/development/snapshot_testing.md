# Snapshot Testing in MultiQC

MultiQC uses [Syrupy](https://github.com/tophat/syrupy) for snapshot testing to ensure that module outputs are documented, reproducible, and consistent. This guide explains how to implement and use snapshot testing for MultiQC modules.

## Overview

Snapshot testing captures the output of your module and stores it as a "snapshot" file. Future test runs compare the current output against the stored snapshot, ensuring that changes to the code don't unexpectedly alter the module's behavior.

### Benefits

1. **Documentation**: Snapshots serve as living documentation of what your module produces
2. **Reproducibility**: Ensures consistent output across different environments
3. **Regression Detection**: Automatically detects when changes affect module output
4. **Table Format Consistency**: Helps maintain consistent table formatting across modules

## Quick Start

### 1. Basic Snapshot Test

The simplest way to add snapshot testing to your module is to inherit from `BaseModuleTest`:

```python
# multiqc/modules/yourmodule/tests/test_yourmodule.py
import pytest
from tests.conftest import BaseModuleTest
from multiqc.modules.yourmodule import MultiqcModule

class TestYourModuleSnapshot(BaseModuleTest):
    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "yourmodule"
    
    # The base class provides these tests automatically:
    # - test_module_data_integrity()
    # - test_module_complete_snapshot()
    # - test_module_raw_data_snapshot()
    # - test_module_general_stats_snapshot()
    # - test_module_sections_snapshot()
```

### 2. Running Tests

```bash
# Run tests and generate initial snapshots
pytest multiqc/modules/yourmodule/tests/ --snapshot-update

# Run tests normally (will compare against snapshots)
pytest multiqc/modules/yourmodule/tests/
```

### 3. Updating Snapshots

When you intentionally change your module's output:

```bash
# Update all snapshots for your module
pytest multiqc/modules/yourmodule/tests/ --snapshot-update

# Update specific test snapshots
pytest multiqc/modules/yourmodule/tests/test_yourmodule.py::TestYourModuleSnapshot::test_module_complete_snapshot --snapshot-update
```

## Advanced Usage

### Custom Snapshot Tests

For more control over what gets snapshotted:

```python
class TestYourModuleCustom:
    def test_specific_parsing_logic(self, data_dir, snapshot, snapshot_config):
        """Test specific parsing logic with custom data."""
        test_file = data_dir / "modules/yourmodule/special_case.log"
        
        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)
        
        # Run module test
        module_snapshot = testing.run_module_test(
            module_class=MultiqcModule,
            data_files=[test_file],
            config_updates=snapshot_config,
        )
        
        # Create custom snapshot
        custom_data = {
            "parsed_metrics": module_snapshot.get_saved_raw_data(),
            "sample_count": len(module_snapshot.get_saved_raw_data().get("multiqc_yourmodule", {})),
            "has_general_stats": bool(module_snapshot.get_general_stats_data()),
        }
        
        assert custom_data == snapshot
```

### Testing Individual Files

Test specific input files separately:

```python
@pytest.mark.parametrize("filename", [
    "sample1.log",
    "sample2.txt",
    "edge_case.log",
])
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
        module_class=MultiqcModule,
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
```

### Testing Table Format Consistency

Ensure your module produces consistent table formats:

```python
def test_table_format_consistency(self, module_snapshot):
    """Test that table headers follow MultiQC standards."""
    headers = module_snapshot.get_general_stats_headers()
    
    for section_key, section_headers in headers.items():
        for header_key, header_config in section_headers.items():
            # Required fields
            assert "title" in header_config, f"Header {header_key} missing title"
            assert "description" in header_config, f"Header {header_key} missing description"
            
            # Title formatting
            title = header_config["title"]
            assert isinstance(title, str), f"Header {header_key} title must be string"
            assert len(title) > 0, f"Header {header_key} title cannot be empty"
            assert title[0].isupper(), f"Header {header_key} title should start with capital letter"
            
            # Description formatting
            description = header_config["description"]
            assert isinstance(description, str), f"Header {header_key} description must be string"
            assert len(description) > 0, f"Header {header_key} description cannot be empty"
```

## Snapshot Data Structure

The `ModuleSnapshot` class captures several types of data:

### Complete Snapshot

```python
{
    "module_info": {
        "name": "Module Name",
        "anchor": "module_anchor",
        "info": "Module description",
        "href": ["https://tool.homepage.com"],
        "doi": ["10.1234/example.doi"]
    },
    "saved_raw_data": {
        "multiqc_modulename": {
            "sample1": {"metric1": 100, "metric2": 200},
            "sample2": {"metric1": 150, "metric2": 250}
        }
    },
    "general_stats_data": {
        "module_anchor": {
            "sample1": {"key_metric": 100},
            "sample2": {"key_metric": 150}
        }
    },
    "general_stats_headers": {
        "module_anchor": {
            "key_metric": {
                "title": "Key Metric",
                "description": "Description of the metric",
                "scale": "YlGn",
                "format": "{:,.0f}"
            }
        }
    },
    "sections": [
        {
            "name": "Section Name",
            "anchor": "section_anchor",
            "description": "Section description",
            "plot_type": "bar_graph",
            "plot_id": "plot_id",
            "plot_data": {...},
            "plot_config": {...}
        }
    ],
    "software_versions": {
        "Tool Name": [("1.2.3", "1.2.3")]
    }
}
```

### Raw Data Only

```python
{
    "multiqc_modulename": {
        "sample1": {"metric1": 100, "metric2": 200},
        "sample2": {"metric1": 150, "metric2": 250}
    }
}
```

### General Stats Only

```python
{
    "data": {
        "module_anchor": {
            "sample1": {"key_metric": 100},
            "sample2": {"key_metric": 150}
        }
    },
    "headers": {
        "module_anchor": {
            "key_metric": {
                "title": "Key Metric",
                "description": "Description of the metric"
            }
        }
    }
}
```

## Best Practices

### 1. Start Simple

Begin with the basic `BaseModuleTest` class and add custom tests as needed:

```python
class TestYourModuleSnapshot(BaseModuleTest):
    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "yourmodule"
```

### 2. Test Edge Cases

Add specific tests for edge cases and unusual input files:

```python
def test_empty_file(self, tmp_path, snapshot, snapshot_config):
    """Test behavior with empty input file."""
    empty_file = tmp_path / "empty.log"
    empty_file.write_text("")
    
    # Test that module handles empty file gracefully
    # (might raise ModuleNoSamplesFound or parse empty data)
```

### 3. Validate Data Integrity

Always include data integrity checks:

```python
def test_data_integrity(self, module_snapshot):
    """Test that parsed data is valid."""
    raw_data = module_snapshot.get_saved_raw_data()
    
    for file_key, file_data in raw_data.items():
        for sample_name, sample_data in file_data.items():
            # Check that all required metrics are present
            required_metrics = ["metric1", "metric2"]
            for metric in required_metrics:
                assert metric in sample_data, f"Missing {metric} for {sample_name}"
            
            # Check that numeric values are valid
            for metric, value in sample_data.items():
                if isinstance(value, (int, float)):
                    assert not (isinstance(value, float) and value != value), f"NaN value for {metric}"
```

### 4. Use Descriptive Test Names

Make test names descriptive so it's clear what failed:

```python
def test_parsing_with_missing_headers(self, snapshot):
    """Test parsing when some expected headers are missing."""
    
def test_parsing_with_malformed_data(self, snapshot):
    """Test parsing when data is malformed or incomplete."""
    
def test_general_stats_table_format(self, snapshot):
    """Test that general stats table follows standard format."""
```

### 5. Group Related Tests

Use classes to group related tests:

```python
class TestYourModuleBasicFunctionality(BaseModuleTest):
    """Basic functionality and snapshot tests."""
    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "yourmodule"

class TestYourModuleEdgeCases:
    """Edge cases and error handling."""
    
class TestYourModuleTableFormat:
    """Table format consistency tests."""
```

## Troubleshooting

### Snapshot Mismatches

When snapshots don't match:

1. **Review the diff**: Syrupy shows exactly what changed
2. **Check if change is intentional**: If yes, update snapshots
3. **Look for environment differences**: Paths, floating-point precision, etc.
4. **Check data normalization**: The `ModuleSnapshot` class normalizes data, but edge cases might slip through

### Cross-Platform Issues

The snapshot system normalizes paths and handles common cross-platform issues, but you might need to:

1. **Exclude problematic keys**: Use `exclude_keys` parameter
2. **Custom normalization**: Override `_normalize_data()` method
3. **Platform-specific snapshots**: Use pytest markers for platform-specific tests

### Large Snapshots

If snapshots become too large:

1. **Test specific components**: Use individual snapshot methods instead of complete snapshots
2. **Filter data**: Create custom snapshots with only relevant data
3. **Split tests**: Break large tests into smaller, focused tests

## Migration Guide

### From Existing Tests

To migrate existing tests to use snapshots:

1. **Keep existing logic tests**: Don't replace all assertions with snapshots
2. **Add snapshot tests**: Use snapshots to complement existing tests
3. **Start with raw data**: Begin by snapshotting `saved_raw_data`
4. **Gradually expand**: Add general stats and sections snapshots

### Example Migration

Before:
```python
def test_module(self, data_dir):
    # Run module
    m = MultiqcModule()
    
    # Check specific values
    assert m.saved_raw_data["multiqc_module"]["sample1"]["metric1"] == 100
    assert len(m.sections) == 2
```

After:
```python
class TestModuleSnapshot(BaseModuleTest):
    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "module"
    
    def test_specific_values(self, module_snapshot):
        """Test specific values (keep existing logic)."""
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data["multiqc_module"]["sample1"]["metric1"] == 100
        assert len(module_snapshot.module.sections) == 2
    
    # Snapshot tests are automatically included from BaseModuleTest
```

## Examples

See these modules for examples of snapshot testing:

- `multiqc/modules/samtools/tests/test_flagstat.py` - Complete example with individual file tests
- `tests/test_module_template.py` - Template with all patterns

## Configuration

### Syrupy Configuration

Syrupy can be configured in `pyproject.toml`:

```toml
[tool.syrupy]
# Include these settings if needed
include_paths = ["tests/", "multiqc/modules/"]
exclude_paths = ["tests/test_module_template.py"]
```

### MultiQC Test Configuration

The `snapshot_config` fixture provides standard configuration:

```python
{
    "preserve_module_raw_data": True,
    "strict": True,
    "make_data_dir": False,
}
```

Override this in your tests if needed:

```python
@pytest.fixture
def custom_snapshot_config(self, snapshot_config):
    return {
        **snapshot_config,
        "yourmodule_custom_option": True,
    }
``` 