---
name: Analyzing Bioinformatics Log Files
description: Examines unfamiliar bioinformatics tool log files to identify structure, data format, parsing strategies, and visualization opportunities. Use when encountering new log formats or designing module parsing logic.
allowed-tools: [Read, Grep, Bash]
---

# Analyzing Bioinformatics Log Files

This skill helps understand and design parsing strategies for bioinformatics tool outputs.

## When to Use

- Examining example files from module requests
- Understanding unfamiliar log file formats
- Designing parsing logic for new modules
- Troubleshooting parsing issues

## Analysis Workflow

### Phase 1: Initial Examination

#### 1.1 File Overview

```bash
# Check file size and type
file example_log.txt
wc -l example_log.txt

# Preview first and last lines
head -20 example_log.txt
tail -20 example_log.txt
```

**Identify**:
- Total lines
- File type (text, JSON, XML, binary)
- Obvious structure patterns
- Section markers

#### 1.2 Format Classification

Determine primary format:

**Structured Formats**:
- JSON: `{` and `}`, key-value pairs
- XML/HTML: `<tags>` and `</tags>`
- YAML: Indented key-value with `:`
- CSV/TSV: Comma or tab-delimited columns

**Semi-Structured Formats**:
- Key-value pairs: `Key: Value`
- INI-style sections: `[Section]`
- Property files: `key=value`

**Unstructured Formats**:
- Free-form text with patterns
- Mixed formats
- Custom delimiters

### Phase 2: Structure Analysis

#### 2.1 Identify Sections

Look for section markers:
```bash
# Common section patterns
grep -n "^##" example_log.txt      # Double hash headers
grep -n "^>" example_log.txt       # Greater-than headers
grep -n "^\[" example_log.txt      # Bracket sections
grep -n "^===" example_log.txt     # Separator lines
```

Document sections found:
```markdown
## File Structure

### Section 1: Header (lines 1-5)
- Format: Key-value pairs
- Contains: Version, timestamp, input file

### Section 2: Statistics (lines 7-15)
- Format: Tab-separated
- Contains: Numeric metrics

### Section 3: Results (lines 17-end)
- Format: Table with headers
- Contains: Per-feature results
```

#### 2.2 Sample Name Detection

Identify where sample name appears:
```bash
# Search for common sample name patterns
grep -i "sample" example_log.txt
grep -i "filename" example_log.txt
grep -i "input" example_log.txt
```

**Common locations**:
- In filename itself (`sample1_tool_output.txt`)
- First line of file
- JSON/YAML key: `"sample": "name"`
- Command line in header

**Extraction strategy**:
```python
# From filename
s_name = filename.replace("_tool_output.txt", "")

# From file content (key-value)
match = re.search(r"Sample:\s*(\S+)", content)
if match:
    s_name = match.group(1)

# From JSON
data = json.loads(content)
s_name = data.get("sample", "unknown")
```

### Phase 3: Metrics Identification

#### 3.1 Locate Key Metrics

Search for numeric values:
```bash
# Find lines with numbers
grep -E "[0-9]+" example_log.txt | head -20

# Find percentages
grep -E "[0-9]+\.[0-9]+%" example_log.txt

# Find scientific notation
grep -E "[0-9]+\.[0-9]+e[+-][0-9]+" example_log.txt
```

#### 3.2 Categorize Metrics

Group by type:

**Count Metrics** (integers):
- Total reads, passed reads, failed reads
- Number of features, variants, alignments
- Best for: Bar charts, general stats

**Rate Metrics** (percentages):
- Pass rate, error rate, coverage
- Best for: General stats (colored by threshold)

**Distribution Metrics** (position → value):
- Quality by position
- Coverage by chromosome
- Best for: Line graphs

**Categorical Metrics** (category → count):
- Read length distribution
- Base composition
- Best for: Bar charts, stacked bar charts

**Continuous Metrics** (x → y):
- GC content distribution
- Insert size distribution
- Best for: Line graphs, density plots

### Phase 4: Parsing Strategy Design

#### 4.1 Choose Parsing Approach

Based on format classification:

**For JSON**:
```python
import json

def parse_json_log(f):
    data = json.loads(f)
    return {
        "metric1": data["stats"]["metric1"],
        "metric2": data["stats"]["metric2"]
    }
```

**For TSV**:
```python
def parse_tsv_log(f):
    data = {}
    for line in f.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            data[parts[0]] = parse_value(parts[1])
    return data
```

**For Key-Value**:
```python
import re

def parse_keyvalue_log(f):
    data = {}
    for line in f.splitlines():
        match = re.match(r"(\S+):\s*(.+)", line)
        if match:
            key = match.group(1).lower().replace(" ", "_")
            value = parse_value(match.group(2))
            data[key] = value
    return data
```

**For Multi-Section**:
```python
def parse_multisection_log(f):
    data = {}
    current_section = None

    for line in f.splitlines():
        # Detect section headers
        if line.startswith("##"):
            current_section = line.strip("#").strip()
            data[current_section] = {}
        # Parse content within section
        elif current_section:
            # Apply section-specific parsing
            pass

    return data
```

**For Mixed Formats**:
```python
def parse_mixed_log(f):
    lines = f.splitlines()
    data = {}

    # Parse header (lines 1-5, key-value)
    for line in lines[:5]:
        match = re.match(r"(\w+):\s*(.+)", line)
        if match:
            data[match.group(1)] = match.group(2)

    # Parse table (lines 7+, TSV)
    for line in lines[7:]:
        parts = line.split("\t")
        # Process table data

    return data
```

#### 4.2 Handle Edge Cases

Consider:
- Empty files
- Missing sections
- Malformed data
- Different versions with different formats
- Optional fields

```python
def parse_robust(f):
    try:
        data = attempt_parsing(f)

        # Validate required fields
        required = ["total_reads", "sample_name"]
        if not all(k in data for k in required):
            log.warning("Missing required fields")
            return None

        return data
    except Exception as e:
        log.error(f"Parsing failed: {e}")
        return None
```

### Phase 5: Visualization Planning

#### 5.1 Match Data to Plot Types

Based on metrics identified:

**Bar Chart** - Best for:
- Category counts (read types, feature classes)
- Pass/fail results
- Comparison across categories

Example:
```python
bargraph.plot(
    {"sample1": {"passed": 1000, "failed": 50}},
    pconfig={"title": "Results Summary"}
)
```

**Line Graph** - Best for:
- Distributions over position
- Quality scores over read length
- Coverage across genome

Example:
```python
linegraph.plot(
    {"sample1": {1: 30, 2: 32, 3: 35}},  # position: value
    pconfig={"title": "Quality by Position", "xlab": "Position", "ylab": "Quality"}
)
```

**Table** - Best for:
- Detailed per-sample statistics
- Many metrics to display
- Sortable data

Example:
```python
table.plot(
    data_by_sample,
    headers={"metric1": {"title": "Metric 1"}}
)
```

**Heatmap** - Best for:
- Matrix data (samples × features)
- Correlation matrices
- Presence/absence across samples

#### 5.2 General Stats Selection

Choose 3-5 most important metrics:
- Represent overall quality/success
- Enable cross-sample comparison
- Include at least one pass/fail metric

**Good choices**:
- Total count (e.g., total reads)
- Success rate (e.g., pass percentage)
- Quality metric (e.g., average quality)
- Key result (e.g., features detected)

**Avoid**:
- Too many metrics (clutters table)
- Redundant metrics (total vs percentage)
- Internal/debug values

### Phase 6: Implementation Checklist

Before implementing module:

- [ ] Format clearly identified
- [ ] Sample name extraction strategy defined
- [ ] Key metrics identified and categorized
- [ ] Parsing approach selected
- [ ] Edge cases considered
- [ ] Visualization types chosen
- [ ] General stats metrics selected
- [ ] Version string location identified (if available)

## Common Patterns Reference

### Pattern 1: FastQC-Style

**Characteristics**:
- Multiple `>>Section` markers
- Mix of table and key-value data
- Pass/Fail indicators

**Example**:
```
>>Basic Statistics  pass
Total Sequences  1000000
Sequence Length  150
%GC  50
>>END_MODULE
```

**Parsing**:
```python
sections = {}
current_section = None

for line in f.splitlines():
    if line.startswith(">>"):
        if "END_MODULE" in line:
            current_section = None
        else:
            current_section = line.strip(">").split()[0]
            sections[current_section] = {}
    elif current_section and "\t" in line:
        key, val = line.split("\t")
        sections[current_section][key] = val
```

### Pattern 2: SAMtools-Style

**Characteristics**:
- Simple TSV or key-value
- One metric per line
- No sections

**Example**:
```
reads mapped: 950000
reads unmapped: 50000
average quality: 35.2
```

**Parsing**:
```python
data = {}
for line in f.splitlines():
    match = re.match(r"(.+?):\s*(.+)", line)
    if match:
        key = match.group(1).replace(" ", "_")
        value = parse_number(match.group(2))
        data[key] = value
```

### Pattern 3: Structured JSON/YAML

**Characteristics**:
- Hierarchical data
- Well-defined schema
- Easy to parse

**Example**:
```json
{
  "sample": "sample1",
  "version": "1.0",
  "metrics": {
    "total": 1000000,
    "passed": 950000
  }
}
```

**Parsing**:
```python
import json
data = json.loads(f)
metrics = data["metrics"]
metrics["sample"] = data["sample"]
return metrics
```

## Related Resources

- [Common Log Patterns](./common-patterns.md)
- [Parsing Examples](./parsing-examples.md)
- [Regular Expression Reference](https://regex101.com/)
- [MultiQC Plot Types](https://multiqc.info/docs/#plots)
