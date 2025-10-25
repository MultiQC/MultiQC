# /new-module Command

This command generates a complete MultiQC module implementation from a GitHub issue containing a module request.

## Overview

When invoked, this command will:
1. Analyze the GitHub issue for module details and example files
2. Generate a complete module implementation following MultiQC patterns
3. Create all necessary configuration updates
4. Generate comprehensive tests
5. Create a pull request with the implementation

## Prerequisites

The target issue must:
- Have the `module: new` label
- Have `priority: high` label (score ≥ 70)
- Include uploaded example files (not copy-pasted text)
- Have complete tool information (name, homepage, description)

## Implementation Process

### Step 1: Issue Analysis

Extract the following information from the issue:
- **Tool Name**: From the "Name of the tool" field
- **Tool Homepage**: Repository URL for popularity/maintenance assessment
- **Tool Description**: Short description for module info
- **Example Files**: Download all uploaded files for analysis
- **Log Pattern**: Any specified filename patterns
- **Expected Plots**: User suggestions for visualizations
- **General Stats**: Suggested metrics for the main table

### Step 2: Example File Analysis

For each uploaded example file:
1. **Download and examine** the file contents
2. **Identify data structure** (tab-separated, JSON, key-value pairs, etc.)
3. **Extract sample data** to understand metrics available
4. **Determine parsing strategy** (regex, split, JSON parsing, etc.)
5. **Identify key metrics** that should appear in general stats
6. **Plan visualizations** based on data types (counts → bar charts, distributions → line plots, etc.)

### Step 3: Module Structure Creation

Create the following directory structure:
```
multiqc/modules/{toolname}/
├── __init__.py
├── {toolname}.py
└── tests/
    ├── __init__.py
    └── test_{toolname}.py
```

### Step 4: Main Module Implementation

Generate `{toolname}.py` following this pattern:

```python
import logging
import re
from typing import Dict, Union, Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table
from multiqc.plots.table_object import TableConfig, ColumnMeta

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    [Tool-specific documentation based on analysis]
    
    [Include any special usage notes, supported versions, etc.]
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="[Tool Name]",
            anchor="[lowercase_tool_name]",
            href="[Tool Homepage URL]",
            info="[Tool Description - must start with capital letter]",
            doi="[DOI if available, otherwise omit this field]",
        )

        # Find and parse log files
        data_by_sample: Dict[str, Dict[str, Union[float, int]]] = {}
        for f in self.find_log_files("[toolname]"):
            parsed_data = self.parse_[toolname]_log(f["f"])
            if parsed_data:
                s_name = f["s_name"]
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                data_by_sample[s_name] = parsed_data
                self.add_data_source(f)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Filter ignored samples
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Add to general stats table
        self.add_general_stats_table(data_by_sample)

        # Add module sections with plots
        self.add_sections(data_by_sample)

        # Write data file (MUST be called at the end!)
        self.write_data_file(data_by_sample, f"multiqc_{toolname}")

    def parse_[toolname]_log(self, log_content: str) -> Dict[str, Union[float, int]]:
        """Parse [toolname] log file and extract metrics."""
        # [Implementation based on example file analysis]
        # Return dictionary of metrics
        pass

    def add_general_stats_table(self, data_by_sample: Dict[str, Dict[str, Union[float, int]]]) -> None:
        """Add key metrics to the general statistics table."""
        headers = {
            # [Define 1-2 visible columns and several hidden ones]
            # Use appropriate color schemes and formats
        }
        self.general_stats_addcols(data_by_sample, headers)

    def add_sections(self, data_by_sample: Dict[str, Dict[str, Union[float, int]]]) -> None:
        """Add module sections with appropriate plots."""
        # [Create sections based on data analysis]
        # Use bargraph.plot() for count data
        # Use linegraph.plot() for distributions/trends
        # Use table.plot() for detailed metrics
        pass
```

### Step 5: Search Pattern Configuration

Add to `multiqc/search_patterns.yaml`:
```yaml
[toolname]:
  fn: "[filename_pattern]"  # Only if tool has standard naming
  contents: "[unique_string_from_logs]"
  num_lines: [appropriate_limit]  # Keep as low as possible for performance
```

### Step 6: Module Registration

Add to `pyproject.toml` in the `[project.entry-points."multiqc.modules.v1"]` section:
```toml
[toolname] = "multiqc.modules.[toolname]:MultiqcModule"
```

### Step 7: Module Ordering

Add to `multiqc/config_defaults.yaml` in the `module_order` list in appropriate position (generally alphabetical unless there's a logical grouping).

### Step 8: Test Implementation

Create `tests/test_{toolname}.py`:
```python
import pytest
from multiqc.modules.[toolname] import MultiqcModule
from multiqc import report, config

def test_[toolname]():
    """Test [toolname] module with provided example files."""
    # [Use actual example files from the issue]
    # Test that module finds files, parses data, and generates expected output
    pass
```

### Step 9: Quality Validation

Ensure the generated module:
- **Follows naming conventions**: lowercase module name, proper class name
- **Has proper error handling**: Raises ModuleNoSamplesFound when appropriate
- **Includes logging**: Uses log.info, log.debug appropriately
- **Manages duplicates**: Handles duplicate sample names
- **Has documentation**: Class docstring with usage notes
- **Follows modern Python**: f-strings, type hints, proper imports
- **Uses double quotes**: For string literals
- **No unnecessary features**: Don't add everything, focus on key metrics

### Step 10: Branch and PR Creation

1. **Create feature branch**: `add-{toolname}-module`
2. **Commit files** with descriptive messages following conventional commits
3. **Create pull request** with:
   - Clear title: "Add [Tool Name] module"
   - Description linking to original issue
   - Summary of implemented features
   - Notes about example file analysis
   - Any limitations or future enhancements

## Code Quality Standards

### Required Patterns
- Always call `self.add_software_version()` even if version not found
- Always call `self.write_data_file()` at the very end
- Always raise `ModuleNoSamplesFound` when no samples found (NOT UserWarning)
- Always use `self.ignore_samples()` to handle sample filtering
- Always add proper entry point in pyproject.toml

### Plotting Guidelines
- Use appropriate plot types for data:
  - **Bar charts**: For counts, categorical data
  - **Line plots**: For distributions, trends over positions
  - **Tables**: For detailed metrics, many columns
  - **Scatter plots**: For correlation data
- Choose meaningful color schemes
- Add proper descriptions and help text
- Make plots interactive when possible

### Performance Considerations
- Set `num_lines` in search patterns as low as possible
- Use efficient parsing (avoid complex regex where possible)  
- Don't load large files unnecessarily
- Consider memory usage with many samples

## Example Analysis Patterns

### Tab-separated Data
```python
def parse_log(self, log_content):
    data = {}
    for line in log_content.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        cols = line.strip().split('\t')
        if len(cols) >= 2:
            data[cols[0]] = float(cols[1])
    return data
```

### Key-Value Format
```python
def parse_log(self, log_content):
    data = {}
    for line in log_content.splitlines():
        if ':' in line:
            key, value = line.split(':', 1)
            try:
                data[key.strip()] = float(value.strip())
            except ValueError:
                continue
    return data
```

### JSON Format
```python
import json

def parse_log(self, log_content):
    try:
        data = json.loads(log_content)
        return {k: v for k, v in data.items() if isinstance(v, (int, float))}
    except json.JSONDecodeError:
        return {}
```

## Success Criteria

A successful module implementation should:
- Parse all provided example files correctly
- Extract meaningful metrics for general stats
- Generate appropriate visualizations
- Follow MultiQC coding standards
- Include comprehensive tests
- Have clear documentation
- Be ready for maintainer review with minimal changes needed

When complete, the module should integrate seamlessly into MultiQC and provide valuable insights to users analyzing the tool's output.