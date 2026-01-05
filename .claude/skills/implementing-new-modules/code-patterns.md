# Code Patterns for MultiQC Modules

## Module Class Structure

### Single-Tool Module

```python
"""MultiQC module to parse output from ToolName"""

import logging
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    ToolName provides [description].

    The module parses output from `toolname command` which produces
    [description of output].
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ToolName",
            anchor="toolname",
            href="https://example.com/toolname",
            info="Brief description of what the tool does.",
            doi="10.xxxx/xxxxx",
        )

        # Parse data
        self.toolname_data = {}
        for f in self.find_log_files("toolname"):
            parsed = self.parse_log(f)
            # ... process parsed data

        # Filter ignored samples
        self.toolname_data = self.ignore_samples(self.toolname_data)

        if len(self.toolname_data) == 0:
            raise ModuleNoSamplesFound

        # Add sections and stats
        self.add_general_stats()
        self.add_sections()

        # Write data file (MUST be at the end)
        self.write_data_file(self.toolname_data, "multiqc_toolname")
```

### Multi-Subtool Module Orchestrator

```python
"""MultiQC module to parse output from ToolName"""

import logging
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .subtool1 import parse_toolname_subtool1
from .subtool2 import parse_toolname_subtool2

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Supported commands:

    - `subtool1`
    - `subtool2`
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ToolName",
            anchor="toolname",
            href="https://example.com/toolname",
            info="Brief description of what the tool does.",
            doi="10.xxxx/xxxxx",
        )

        n = dict()

        n["subtool1"] = parse_toolname_subtool1(self)
        if n["subtool1"] > 0:
            log.info(f"Found {n['subtool1']} subtool1 reports")

        n["subtool2"] = parse_toolname_subtool2(self)
        if n["subtool2"] > 0:
            log.info(f"Found {n['subtool2']} subtool2 reports")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
```

### Submodule Parser Function

```python
"""MultiQC submodule to parse output from toolname subtool"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule, config
from multiqc.plots import table, bargraph

log = logging.getLogger(__name__)


def parse_toolname_subtool(module: BaseMultiqcModule) -> int:
    """Find toolname subtool logs and parse their data"""

    data: Dict[str, Dict] = {}

    for f in module.find_log_files("toolname/subtool"):
        parsed = parse_report()
        for s_name, sample_data in parsed.items():
            s_name = module.clean_s_name(s_name, f)  # Only needed if not using f['s_name']
            if s_name in data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="subtool")
            data[s_name] = sample_data

            # Required call - even if version is None
            module.add_software_version(None, s_name)

    # Filter ignored samples
    data = module.ignore_samples(data)

    if len(data) == 0:
        return 0

    # Add general stats
    add_general_stats(module, data)

    # Add detailed section
    add_section(module, data)

    # Write data file
    module.write_data_file(data, "multiqc_toolname_subtool")

    return len(data)
```

## Parsing Patterns

### Key-Value Pair Parsing

```python
def parse_key_value_report(file_content: str) -> Dict:
    """Parse key: value format output."""
    data = {}
    for line in file_content.strip().split("\n"):
        if ":" in line:
            key, value = line.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
            try:
                data[key] = float(value.replace(",", "").rstrip("%"))
            except ValueError:
                data[key] = value
    return data
```

### JSON Parsing

```python
import json

def parse_json_report(file_content: str) -> Dict:
    """Parse JSON format output."""
    try:
        return json.loads(file_content)
    except json.JSONDecodeError:
        return {}
```

## General Stats Patterns

### Standard Headers Definition

```python
general_stats_headers = {
    "read_count": {
        "title": "# Reads",
        "description": f"Total read count ({config.read_count_desc})",
        "scale": "Blues",
        "shared_key": "read_count",
    },
    "base_count": {
        "title": "Total bp",
        "description": f"Total bases ({config.base_count_desc})",
        "scale": "Greens",
        "shared_key": "base_count",
        "hidden": True,
    },
    "percentage_metric": {
        "title": "Metric%",
        "description": "Description of metric",
        "min": 0,
        "max": 100,
        "scale": "RdYlGn",  # Use for quality (higher=better)
        "suffix": "%",
    },
    "gc_content": {
        "title": "GC%",
        "description": "GC content percentage",
        "min": 0,
        "max": 100,
        "scale": "RdYlBu",  # Use for middle-is-best metrics
        "suffix": "%",
    },
}

# Get headers with config integration
headers = module.get_general_stats_headers(all_headers=general_stats_headers)

# Add to general stats table
if headers:
    module.general_stats_addcols(data, headers, namespace="toolname")
```

## Visualization Patterns

### Table Plot

```python
from multiqc.plots import table

table_headers = {
    "column1": {
        "title": "Column 1",
        "description": "Description",
        "format": "{:,.0f}", # Only needed for integers
        "scale": "Blues",
    },
    # ... more columns
    # Headers without data, or data without headers, will be ignored - no conditionals required
}

module.add_section(
    name="Detailed Stats",
    anchor="toolname-stats",
    description="Detailed statistics from <code>toolname</code>.",
    plot=table.plot(
        data,
        table_headers,
        pconfig={
            "id": "toolname-stats-table",
            "title": "ToolName: Statistics",
            "namespace": "toolname",
        },
    ),
)
```

### Bar Graph

```python
from multiqc.plots import bargraph

# Only needed if data needs to be reshaped.
bargraph_data = {
    s_name: {
        "Category1": d.get("count1", 0),
        "Category2": d.get("count2", 0),
    }
    for s_name, d in data.items()
}

module.add_section(
    name="Counts",
    anchor="toolname-counts",
    description="Count data from toolname.",
    plot=bargraph.plot(
        bargraph_data,
        pconfig={
            "id": "toolname-counts-plot",
            "title": "ToolName: Counts",
            "ylab": "Count",
            "cpswitch": False,  # Disable count/percentage switch
        },
    ),
)
```

### Line Graph

```python
from multiqc.plots import linegraph

# Data format: {sample: {x1: y1, x2: y2, ...}}
line_data = {
    s_name: {pos: val for pos, val in enumerate(d.get("values", []))}
    for s_name, d in data.items()
}

module.add_section(
    name="Distribution",
    anchor="toolname-dist",
    description="Value distribution.",
    plot=linegraph.plot(
        line_data,
        pconfig={
            "id": "toolname-dist-plot",
            "title": "ToolName: Distribution",
            "xlab": "Position",
            "ylab": "Value",
        },
    ),
)
```

## **init**.py Pattern

```python
from .toolname import MultiqcModule

__all__ = ["MultiqcModule"]
```
