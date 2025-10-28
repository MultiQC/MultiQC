---
name: Creating MultiQC Modules
description: Generates complete MultiQC module implementations from GitHub issue specifications, including parsing code, visualizations, tests, and configuration. Use when creating new bioinformatics tool modules or responding to module requests.
allowed-tools: [Read, Write, Edit, Bash, Glob, Grep, WebFetch]
---

# Creating MultiQC Modules

Generate production-ready MultiQC modules following project standards.

## When to Use

- Creating new modules from GitHub issues
- Converting example log files into module code
- Implementing module requests

## Prerequisites

Verify before starting:

- [ ] Issue has `module: new` label
- [ ] Example files uploaded (not pasted text)
- [ ] Tool name and homepage provided
- [ ] Priority score ≥ 70 (or maintainer approval)

## Critical Requirements

**Must follow** (enforced by linting):

- ✅ Raise `ModuleNoSamplesFound` when no samples found
- ✅ Call `self.add_software_version()` even if version unknown
- ✅ Call `self.write_data_file()` at the END
- ✅ Module `info` starts with capital letter
- ✅ Use double quotes for strings
- ✅ No shebang lines
- ✅ Modern Python 3 syntax (f-strings, no `__future__`)

## Implementation Workflow

### Phase 1: Analysis

**Extract from issue**:

- Tool name (lowercase-hyphenated for code, Display Name for UI)
- Homepage URL (for `href`)
- One-line description (for `info`, starts with capital)
- Example file URLs

**Analyze example files** using the log analysis skill:

- File format and structure
- Sample name location
- Key metrics to extract
- Suggested parsing strategy

### Phase 2: Generate Module

#### Directory Structure

```bash
mkdir -p multiqc/modules/{toolname}/tests/data
touch multiqc/modules/{toolname}/__init__.py
touch multiqc/modules/{toolname}/{toolname}.py
touch multiqc/modules/{toolname}/tests/__init__.py
touch multiqc/modules/{toolname}/tests/test_{toolname}.py
```

#### Main Module File

**Use the template**: [templates/module_template.py](./templates/module_template.py)

Replace placeholders:

- `{TOOLNAME}` / `{toolname}` → actual tool name
- `{Display Name}` → human-readable name
- `{homepage_url}` → repository/docs URL
- `{One-line description...}` → from issue

**Key implementation points**:

1. Parse files in `__init__` using `find_log_files()`
2. Filter with `ignore_samples()`
3. Always call `add_software_version()` (even with None)
4. Add general stats for 3-5 key metrics
5. Create appropriate visualizations
6. Call `write_data_file()` LAST

See [examples/simple-module.md](./examples/simple-module.md) or [examples/complex-module.md](./examples/complex-module.md) for patterns.

#### Configuration Updates

**`multiqc/search_patterns.yaml`**:

```yaml
{ toolname }:
  fn: "*_{toolname}.txt" # Adjust pattern to match actual files
```

**`pyproject.toml`** entry point:

```toml
[project.entry-points."multiqc.modules.v1"]
{toolname} = "multiqc.modules.{toolname}.{toolname}"
```

### Phase 3: Testing

**Add test data**:
Copy example files to `multiqc/modules/{toolname}/tests/data/`

**Write test** in `tests/test_{toolname}.py`:

```python
from multiqc.modules.{toolname} import MultiqcModule

def test_{toolname}(data_dir):
    module = MultiqcModule()
    assert len(module.data_by_sample) > 0
    # Verify key metrics present
```

**Run validation**:

```bash
pytest multiqc/modules/{toolname}/tests/ -v
pre-commit run --all-files
```

### Phase 4: Submit

**Create branch and commit**:

```bash
git checkout -b {toolname}-module
git add multiqc/modules/{toolname}/ multiqc/search_patterns.yaml pyproject.toml
git commit -m "feat: add {Tool Name} module

Implements MultiQC module for {Tool Name} (#{issue_number}).

- Parse {command} output
- Extract key metrics
- Generate visualizations
- Full test coverage

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

**Create PR**:

```bash
git push -u origin {toolname}-module
gh pr create --title "Add {Tool Name} module" --body "
## Summary
Implements MultiQC module for {Tool Name} (#{issue_number}).

### Features
- Parses output files
- Key metrics in general stats
- Visualizations for data distribution
- Comprehensive tests

Closes #{issue_number}

🤖 Generated with [Claude Code](https://claude.ai/code)"
```

## Quality Checklist

**Before submitting**:

### Standards Compliance

- [ ] Raises `ModuleNoSamplesFound` when appropriate
- [ ] Calls `add_software_version()`
- [ ] Calls `write_data_file()` last
- [ ] Entry point in `pyproject.toml`
- [ ] Search pattern in `search_patterns.yaml`

### Code Quality

- [ ] No linting errors (`ruff check`)
- [ ] Type hints pass (`mypy`)
- [ ] Double quotes, f-strings
- [ ] No unused imports

### Functionality

- [ ] Parses all example files correctly
- [ ] Handles malformed data gracefully
- [ ] Sample names extracted properly
- [ ] Appropriate metrics in general stats (3-5 key ones)

### Testing

- [ ] Tests pass (`pytest`)
- [ ] Coverage >80%
- [ ] Uses real example data
- [ ] Tests edge cases

### Documentation

- [ ] Module docstring complete
- [ ] Plot titles clear
- [ ] General stats metrics described

## Plot Selection Guide

**Choose based on data type**:

- **Bar Graph**: Categories, pass/fail, counts by group

  ```python
  bargraph.plot(data, pconfig={"title": "..."})
  ```

- **Line Graph**: Distributions, quality by position, trends

  ```python
  linegraph.plot(data, pconfig={"xlab": "Position", "ylab": "Quality"})
  ```

- **Table**: Detailed per-sample metrics, many values
  ```python
  table.plot(data, headers)
  ```

See examples for complete implementations.

## General Stats Selection

**Choose 3-5 metrics that**:

- Represent overall quality/success
- Enable cross-sample comparison
- Include at least one pass/fail metric

**Good choices**:

- Total count (reads, features, variants)
- Success rate (percentage passed)
- Quality metric (average score)
- Key result (features detected)

**Avoid**:

- Too many metrics (>5)
- Redundant metrics (total + percentage)
- Internal/debug values

## Troubleshooting

**Module not found**:

- Check entry point in `pyproject.toml`
- Reinstall: `pip install -e .`

**Files not detected**:

- Verify pattern in `search_patterns.yaml`
- Test with actual file: `multiqc . --module {toolname}`

**Tests failing**:

- Check example data in `tests/data/`
- Verify parsing handles actual file format
- Run with `-vv` for details

## Related Resources

- [Module Template](./templates/module_template.py) - Complete code structure
- [Simple Module Example](./examples/simple-module.md) - TSV parsing
- [Complex Module Example](./examples/complex-module.md) - JSON with multiple plots
- [Log Analysis Skill](../analyzing-log-files/SKILL.md) - Format identification
- [MultiQC Guidelines](../../CLAUDE.md) - Project standards
- [Plot Types](https://multiqc.info/docs/#plots) - Official documentation
