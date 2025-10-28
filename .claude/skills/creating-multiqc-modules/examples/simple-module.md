# Example: Simple MultiQC Module

This example shows a straightforward module implementation for a tool with simple TSV output.

## Scenario

**Tool**: `readstats` - A simple tool that counts reads in FASTQ files

**Output Format** (`readstats_output.txt`):

```
Total Reads: 1000000
Passed QC: 950000
Failed QC: 50000
Average Length: 150
```

## Implementation

### 1. Module Code

```python
import logging
import re
from typing import Dict, Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC module to parse output from **ReadStats**.

    ReadStats counts reads and performs basic quality filtering on FASTQ files.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ReadStats",
            anchor="readstats",
            href="https://github.com/example/readstats",
            info="Counts reads and performs quality filtering",
        )

        data_by_sample = {}

        for f in self.find_log_files("readstats"):
            parsed_data = self.parse_readstats_log(f["f"])
            if parsed_data:
                s_name = f["s_name"]
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found: {s_name}")
                data_by_sample[s_name] = parsed_data
                self.add_data_source(f)

        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            raise ModuleNoSamplesFound

        self.add_software_version(None)
        log.info(f"Found {len(data_by_sample)} reports")

        # Add to general stats
        headers = {
            "total_reads": {
                "title": "Total Reads",
                "description": "Total number of reads in FASTQ file",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "passed_qc": {
                "title": "Passed QC",
                "description": "Number of reads passing quality filters",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "pass_rate": {
                "title": "Pass Rate",
                "description": "Percentage of reads passing QC",
                "scale": "RdYlGn",
                "format": "{:,.1f}%",
                "max": 100,
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

        # Write data file
        self.write_data_file(data_by_sample, "multiqc_readstats")

    def parse_readstats_log(self, f: str) -> Dict[str, Union[int, float]]:
        """Parse readstats output file."""
        data = {}

        for line in f.splitlines():
            # Match pattern "Key: Value"
            match = re.match(r"(.+?):\s*(\d+(?:\.\d+)?)", line)
            if match:
                key = match.group(1).lower().replace(" ", "_")
                value = float(match.group(2))
                data[key] = value

        # Calculate pass rate
        if "total_reads" in data and "passed_qc" in data:
            data["pass_rate"] = (data["passed_qc"] / data["total_reads"]) * 100

        return data if data else None
```

### 2. Search Pattern

Add to `multiqc/search_patterns.yaml`:

```yaml
readstats:
  fn: "*readstats_output.txt"
```

### 3. Entry Point

Add to `pyproject.toml`:

```toml
[project.entry-points."multiqc.modules.v1"]
readstats = "multiqc.modules.readstats.readstats"
```

### 4. Test

```python
from multiqc.modules.readstats import MultiqcModule


def test_readstats():
    """Test readstats module."""
    module = MultiqcModule()

    assert len(module.data_by_sample) == 1
    assert "sample1" in module.data_by_sample

    data = module.data_by_sample["sample1"]
    assert data["total_reads"] == 1000000
    assert data["passed_qc"] == 950000
    assert data["pass_rate"] == 95.0
```

## Key Takeaways

1. **Simple parsing**: Use regex for key-value format
2. **Calculated metrics**: Derive pass_rate from raw counts
3. **No plots needed**: Simple stats work well in general stats table
4. **Clean implementation**: ~80 lines of code
