# Module Structure Templates

## Single-Tool Module

### Directory Structure

```
multiqc/modules/toolname/
├── __init__.py
├── toolname.py
└── tests/
    ├── __init__.py
    └── test_toolname.py
```

### **init**.py

```python
from .toolname import MultiqcModule

__all__ = ["MultiqcModule"]
```

### toolname.py Template

```python
"""MultiQC module to parse output from ToolName"""

import logging
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table, bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    ToolName is a tool for [description].

    The module parses [output type] from `toolname [command]`.

    Supported output formats:
    - Format 1
    - Format 2 (with `--option`)
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ToolName",
            anchor="toolname",
            href="https://example.com/toolname",
            info="Brief description starting with capital letter.",
            doi="10.xxxx/journal.xxxxx",
        )

        self.toolname_data: Dict[str, Dict] = {}

        for f in self.find_log_files("toolname"):
            parsed = self._parse_log(f)
            for s_name, data in parsed.items():
                s_name = self.clean_s_name(s_name, f) # If needed
                if s_name in self.toolname_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name=s_name)
                self.toolname_data[s_name] = data
                self.add_software_version(None, s_name)

        self.toolname_data = self.ignore_samples(self.toolname_data)

        if len(self.toolname_data) == 0:
            raise ModuleNoSamplesFound

        self._add_general_stats()
        self._add_stats_table()

        self.write_data_file(self.toolname_data, "multiqc_toolname")

    def _parse_log(self, f):
        """Parse toolname output file."""
        # Implementation here
        pass

    def _add_general_stats(self):
        """Add columns to general stats table."""
        headers = {
            "metric1": {
                "title": "Metric 1",
                "description": "Description",
                "scale": "Blues",
            },
        }
        headers = self.get_general_stats_headers(all_headers=headers)
        if headers:
            self.general_stats_addcols(self.toolname_data, headers)

    def _add_stats_table(self):
        """Add detailed stats table section."""
        headers = {
            "metric1": {"title": "Metric 1", "description": "Description"},
        }
        self.add_section(
            name="Statistics",
            anchor="toolname-stats",
            description="Statistics from <code>toolname</code>.",
            plot=table.plot(
                self.toolname_data,
                headers,
                pconfig={"id": "toolname-table", "title": "ToolName: Stats"},
            ),
        )
```

---

## Multi-Subtool Module

### Directory Structure

```
multiqc/modules/toolname/
├── __init__.py
├── toolname.py          # Orchestrator
├── subtool1.py          # First subtool parser
├── subtool2.py          # Second subtool parser
└── tests/
    ├── __init__.py
    ├── test_subtool1.py
    └── test_subtool2.py
```

### **init**.py

```python
from .toolname import MultiqcModule

__all__ = ["MultiqcModule"]
```

### toolname.py (Orchestrator)

```python
"""MultiQC module to parse output from ToolName"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .subtool1 import parse_toolname_subtool1
from .subtool2 import parse_toolname_subtool2

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    ToolName is a toolkit for [description].

    Supported commands:

    - `subtool1`
    - `subtool2`

    #### subtool1

    Description of subtool1 and how to generate output.

    #### subtool2

    Description of subtool2 and how to generate output.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ToolName",
            anchor="toolname",
            href="https://example.com/toolname",
            info="Toolkit for [description].",
            doi="10.xxxx/journal.xxxxx",
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

### subtool1.py (Submodule)

```python
"""MultiQC submodule to parse output from toolname subtool1"""

import logging
from typing import Dict, Optional

from multiqc import BaseMultiqcModule, config
from multiqc.plots import table, bargraph

log = logging.getLogger(__name__)


def parse_toolname_subtool1(module: BaseMultiqcModule) -> int:
    """Find toolname subtool1 logs and parse their data"""

    data: Dict[str, Dict] = {}

    for f in module.find_log_files("toolname/subtool1"):
        parsed = parse_report(f["f"], f["s_name"])
        for s_name, sample_data in parsed.items():
            s_name = module.clean_s_name(s_name, f)
            if s_name in data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="subtool1")
            data[s_name] = sample_data
            module.add_software_version(None, s_name)

    data = module.ignore_samples(data)

    if len(data) == 0:
        return 0

    # Add general stats
    headers = get_general_stats_headers()
    stats_headers = module.get_general_stats_headers(all_headers=headers)
    if stats_headers:
        module.general_stats_addcols(data, stats_headers, namespace="subtool1")

    # Add section
    add_section(module, data)

    # Write data
    module.write_data_file(data, "multiqc_toolname_subtool1")

    return len(data)


def parse_report(file_content: str, fallback_name: Optional[str] = None) -> Dict[str, Dict]:
    """Parse subtool1 output file."""
    parsed_data: Dict[str, Dict] = {}
    # Implementation here
    return parsed_data


def get_general_stats_headers() -> Dict:
    """Return general stats header definitions."""
    return {
        "metric1": {
            "title": "Metric 1",
            "description": "Description",
            "scale": "Blues",
        },
    }


def add_section(module: BaseMultiqcModule, data: Dict[str, Dict]):
    """Add detailed section to report."""
    headers = {
        "metric1": {"title": "Metric 1", "description": "Description"},
    }
    module.add_section(
        name="Subtool1",
        anchor="toolname-subtool1",
        description="Output from <code>toolname subtool1</code>.",
        plot=table.plot(
            data,
            headers,
            pconfig={
                "id": "toolname-subtool1-table",
                "title": "ToolName: Subtool1",
                "namespace": "toolname",
            },
        ),
    )
```

---

## Test File Template

### tests/test_subtool.py

```python
"""Tests for the toolname subtool module"""

import pytest

from multiqc.modules.toolname.subtool import parse_report


# Sample output with all columns
SAMPLE_FULL = """header1\theader2\theader3
value1\tvalue2\tvalue3
"""

# Sample with minimal columns
SAMPLE_MINIMAL = """header1\theader2\theader3
value1\tvalue2\tvalue3
"""

# Empty file
SAMPLE_EMPTY = """header1\theader2\theader3
"""

# Invalid format
SAMPLE_INVALID = """wrong\theaders
data\there
"""


class TestParseReport:
    """Tests for parse_report function"""

    def test_parse_full_output(self):
        """Test parsing output with all columns"""
        result = parse_report(SAMPLE_FULL)
        assert len(result) == 1
        # Add specific assertions

    def test_parse_minimal_output(self):
        """Test parsing minimal output"""
        result = parse_report(SAMPLE_MINIMAL)
        assert len(result) == 1

    def test_parse_empty_file(self):
        """Test parsing empty file returns empty dict"""
        result = parse_report(SAMPLE_EMPTY)
        assert result == {}

    def test_parse_invalid_format(self):
        """Test parsing invalid format returns empty dict"""
        result = parse_report(SAMPLE_INVALID)
        assert result == {}

    def test_fallback_sample_name(self):
        """Test fallback sample name for stdin input"""
        stdin_data = """header1\theader2\theader3
-\tvalue2\tvalue3
"""
        result = parse_report(stdin_data, fallback_sample_name="my_sample")
        assert "my_sample" in result
```

---

## Registration Files

### search_patterns.yaml Entry

```yaml
# For single-tool module
toolname:
  contents: "Expected header or unique string"
  num_lines: 1

# For multi-subtool module
toolname/subtool1:
  contents_re: "^header1\\s+header2\\s+header3"
  num_lines: 1

toolname/subtool2:
  fn: "*subtool2*.txt"
```

### pyproject.toml Entry

```toml
[project.entry-points."multiqc.modules.v1"]
# ... other modules ...
toolname = "multiqc.modules.toolname:MultiqcModule"
# ... other modules ...
```

Note: Entry points must be in alphabetical order.
