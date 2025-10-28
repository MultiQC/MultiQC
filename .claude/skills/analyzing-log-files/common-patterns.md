# Common Bioinformatics Log File Patterns

Quick reference for frequently encountered log file formats in bioinformatics tools.

## Format Categories

### 1. Simple Key-Value

**Pattern**: `Key: Value` or `Key = Value`

**Example**:

```
Total reads: 1000000
Mapped reads: 950000
Average quality: 35.2
```

**Regex**: `r"(.+?)[:=]\s*(.+)"`

**Python Parsing**:

```python
data = {}
for line in f.splitlines():
    if ":" in line or "=" in line:
        sep = ":" if ":" in line else "="
        parts = line.split(sep, 1)
        key = parts[0].strip().lower().replace(" ", "_")
        value = parts[1].strip()
        data[key] = parse_number(value)
```

---

### 2. Tab-Separated Values (TSV)

**Pattern**: Columns separated by tabs

**Example**:

```
Feature	Count	Percentage
Mapped	950000	95.0
Unmapped	50000	5.0
```

**Python Parsing**:

```python
data = {}
lines = f.splitlines()
headers = lines[0].split("\t")

for line in lines[1:]:
    parts = line.split("\t")
    row_data = dict(zip(headers, parts))
    # Process row
```

---

### 3. Multi-Section Format

**Pattern**: Sections marked by headers (`##`, `>`, `[Section]`)

**Example**:

```
##Basic Statistics
Total Sequences: 1000000
Sequence Length: 150

##Quality Scores
Position	Mean	Median
1	35	36
2	34	35
```

**Python Parsing**:

```python
sections = {}
current_section = None

for line in f.splitlines():
    if line.startswith("##"):
        current_section = line.strip("#").strip()
        sections[current_section] = {}
    elif current_section and ":" in line:
        key, val = line.split(":", 1)
        sections[current_section][key.strip()] = val.strip()
```

---

### 4. JSON Format

**Pattern**: Hierarchical JSON structure

**Example**:

```json
{
  "sample": "sample1",
  "version": "1.0.0",
  "stats": {
    "total_reads": 1000000,
    "mapped_reads": 950000
  }
}
```

**Python Parsing**:

```python
import json

data = json.loads(f)
stats = data["stats"]
stats["sample"] = data["sample"]
stats["version"] = data.get("version")
return stats
```

---

### 5. YAML Format

**Pattern**: Indented key-value pairs

**Example**:

```yaml
sample: sample1
version: 1.0.0
stats:
  total_reads: 1000000
  mapped_reads: 950000
```

**Python Parsing**:

```python
import yaml

data = yaml.safe_load(f)
stats = data["stats"]
stats["sample"] = data["sample"]
return stats
```

---

### 6. Space-Separated Columns

**Pattern**: Columns separated by spaces

**Example**:

```
Chr  Start    End      Coverage
chr1 1000     2000     45.2
chr1 2000     3000     52.1
```

**Python Parsing**:

```python
data = []
lines = f.splitlines()
headers = lines[0].split()

for line in lines[1:]:
    parts = line.split()
    row = dict(zip(headers, parts))
    data.append(row)
```

---

### 7. Log-Style Format

**Pattern**: Timestamp/level followed by message

**Example**:

```
[2024-01-15 10:30:15] INFO: Processing sample1
[2024-01-15 10:30:16] INFO: Total reads: 1000000
[2024-01-15 10:30:17] INFO: Mapped reads: 950000
```

**Python Parsing**:

```python
import re

data = {}
for line in f.splitlines():
    match = re.search(r"(\w[\w\s]+):\s*(\d+)", line)
    if match:
        key = match.group(1).lower().replace(" ", "_")
        value = int(match.group(2))
        data[key] = value
```

---

### 8. FASTA-Style Headers

**Pattern**: `>header` followed by data

**Example**:

```
>sample1
Total=1000000
Mapped=950000
>sample2
Total=800000
Mapped=750000
```

**Python Parsing**:

```python
data = {}
current_sample = None

for line in f.splitlines():
    if line.startswith(">"):
        current_sample = line[1:].strip()
        data[current_sample] = {}
    elif current_sample and "=" in line:
        key, val = line.split("=")
        data[current_sample][key] = parse_number(val)
```

---

### 9. XML Format

**Pattern**: XML/HTML tags

**Example**:

```xml
<sample name="sample1">
  <stats>
    <total_reads>1000000</total_reads>
    <mapped_reads>950000</mapped_reads>
  </stats>
</sample>
```

**Python Parsing**:

```python
import xml.etree.ElementTree as ET

root = ET.fromstring(f)
sample_name = root.get("name")
stats = {}

for stat in root.find("stats"):
    stats[stat.tag] = parse_number(stat.text)

stats["sample"] = sample_name
return stats
```

---

### 10. Fixed-Width Columns

**Pattern**: Columns at fixed positions

**Example**:

```
Feature     Count      Percentage
Mapped      950000     95.0
Unmapped    50000      5.0
```

**Python Parsing**:

```python
# Define column positions
columns = [
    ("feature", 0, 12),
    ("count", 12, 23),
    ("percentage", 23, 35)
]

data = []
for line in lines[1:]:  # Skip header
    row = {}
    for name, start, end in columns:
        value = line[start:end].strip()
        row[name] = parse_number(value)
    data.append(row)
```

## Special Patterns

### Version Extraction

Common locations:

```python
# First line
version = lines[0].split()[-1]

# Comment header
match = re.search(r"version[:\s]+(\S+)", f, re.IGNORECASE)
if match:
    version = match.group(1)

# JSON/YAML field
data = json.loads(f)
version = data.get("version")
```

### Sample Name Extraction

```python
# From filename
s_name = filename.split("_")[0]
s_name = filename.replace("_output.txt", "")

# From file content
match = re.search(r"[Ss]ample[:\s]+(\S+)", f)
if match:
    s_name = match.group(1)

# From JSON
data = json.loads(f)
s_name = data.get("sample", data.get("id", "unknown"))
```

### Number Parsing

```python
def parse_number(value):
    """Convert string to int or float."""
    # Remove common units and formatting
    value = value.replace(",", "").strip()
    value = re.sub(r"[%]", "", value)

    try:
        # Try int first
        return int(value)
    except ValueError:
        try:
            # Then float
            return float(value)
        except ValueError:
            # Return as string if not a number
            return value
```

## Detection Tips

### Quick Format Detection

```python
def detect_format(f):
    """Detect log file format."""
    first_line = f.splitlines()[0] if f else ""

    if first_line.startswith("{"):
        return "json"
    elif first_line.startswith("<"):
        return "xml"
    elif "\t" in f and f.count("\t") > 10:
        return "tsv"
    elif ":" in f and f.count(":") > 5:
        return "key-value"
    elif f.startswith("##") or ">>END_MODULE" in f:
        return "multi-section"
    else:
        return "unknown"
```

### Structure Analysis

```bash
# Count common delimiters
grep -o ":" file.txt | wc -l    # Colons
grep -o "\t" file.txt | wc -l   # Tabs
grep -o "," file.txt | wc -l    # Commas

# Find section markers
grep -E "^(##|>>|\[)" file.txt

# Look for JSON/XML
head -1 file.txt  # Check first character
```

## Best Practices

1. **Always validate input**: Check for required fields
2. **Handle malformed data**: Use try-except blocks
3. **Log warnings**: Help users debug issues
4. **Return None on failure**: Don't raise exceptions in parsing
5. **Test with multiple files**: Ensure robustness
6. **Document assumptions**: Comment expected format
