# Implementation Checklist

## Phase 1: Research and Planning

### Understand the Tool

- [ ] Read tool documentation for output format specification
- [ ] Identify all output files the tool produces
- [ ] Determine if tool has multiple subcommands with different outputs
- [ ] Note any version-specific output format differences

### Gather Test Data

- [ ] Check `MultiQC/test-data/data/modules/` for existing test files
- [ ] If not present, obtain sample output files
- [ ] Ensure test data covers:
  - [ ] Standard output format
  - [ ] Output with all optional columns (e.g., `--all` flag)
  - [ ] Edge cases (empty, minimal data)

### Choose Architecture

- [ ] Single-tool module (one output format)
- [ ] Multi-subtool module (multiple subcommands)
- [ ] Reference similar existing module for patterns

## Phase 2: Module Structure

### Create Directory Structure

```bash
mkdir -p multiqc/modules/toolname/tests
touch multiqc/modules/toolname/__init__.py
touch multiqc/modules/toolname/toolname.py
touch multiqc/modules/toolname/tests/__init__.py
```

### For Multi-subtool Modules

```bash
touch multiqc/modules/toolname/subtool.py
touch multiqc/modules/toolname/tests/test_subtool.py
```

## Phase 3: Parser Implementation

### Main Module Class (`toolname.py`)

- [ ] Import `BaseMultiqcModule` and `ModuleNoSamplesFound`
- [ ] Create `MultiqcModule` class with docstring describing supported commands
- [ ] Call `super().__init__()` with:
  - [ ] `name` - Display name
  - [ ] `anchor` - URL anchor (lowercase, no spaces)
  - [ ] `href` - Tool homepage URL
  - [ ] `info` - Brief description (starts with capital letter)
  - [ ] `doi` - Publication DOI if available

### Submodule Parser (for multi-subtool)

- [ ] Create `parse_toolname_subtool(module)` function
- [ ] Return count of samples found (not 0 if successful)

### Core Parser Logic

- [ ] Use `module.find_log_files("toolname/subtool")` to discover files
- [ ] Handle both tab-separated and space-separated formats if applicable
- [ ] Use `dict(zip(headers, values))` for clean parsing
- [ ] Convert numeric values appropriately (int vs float)
- [ ] Rename problematic column names (e.g., `Q20(%)` → `Q20_pct`)

### Sample Name Handling

- [ ] Use `module.clean_s_name(s_name, f)` for cleanup if not using `f["s_name"]` (ie. coming from file contents)
- [ ] Add to `fn_clean_exts` in `config_defaults.yaml` if tool has a standard extension that needs removing

### Required Calls

- [ ] `module.add_data_source(f, s_name=s_name, section="subtool")`
- [ ] `module.add_software_version(version, s_name)` (even if version is None)
- [ ] `module.ignore_samples(data_dict)` to filter ignored samples
- [ ] `module.write_data_file(data, "multiqc_toolname_subtool")` at the END

### Error Handling

- [ ] Raise `ModuleNoSamplesFound` when no samples found (NOT `UserWarning`)
- [ ] Use `log.debug()` for skipped files
- [ ] Handle malformed input gracefully, with `log.debug()` messages

## Phase 4: Visualizations

### General Stats Table

- [ ] Define headers dict with appropriate keys
- [ ] Ensure number of headers is apprpriate (not too many)
- [ ] Include for each column:
  - [ ] `title` - Short display title
  - [ ] `description` - Tooltip text (use `config.read_count_desc` etc. where appropriate)
  - [ ] `scale` - Color scale name
  - [ ] `shared_key` - For related columns (e.g., `"read_count"`, `"base_count"`)
  - [ ] `hidden` - True for less important columns
  - [ ] `min`/`max` - For percentage columns
  - [ ] `suffix` - Units (e.g., `"%"`, `" bp"`)
- [ ] Check that only keys with non-default values are included
- [ ] Use `module.get_general_stats_headers()` for config integration
- [ ] Call `module.general_stats_addcols(data, headers, namespace="toolname")`

### Detailed Section

- [ ] Add table with `table.plot()` for detailed metrics
- [ ] Add bar graphs with `bargraph.plot()` for counts/lengths
- [ ] Add line plots if applicable, or heatmaps, violin, box, scatter plots
- [ ] Use `module.add_section()` with:
  - [ ] `name` - Section title
  - [ ] `anchor` - Unique anchor ID
  - [ ] `description` - Section description
  - [ ] `plot` - Plot object

### Color Scale Guidelines

Use scales from ColorBrewer:

- [ ] Sequential
  - [ ] OrRd, PuBu, BuPu, Oranges, BuGn, YlOrBr, YlGn, Reds, RdPu, Greens, YlGnBu, Purples, GnBu, Greys, YlOrRd, PuRd, Blues, PuBuGn
- [ ] Diverging
  - [ ] Spectral, RdYlGn, RdBu, PiYG, PRGn, RdYlBu, BrBG, RdGy, PuOr
- [ ] Qualitative
  - [ ] Set2, Accent, Set1, Set3, Dark2, Paired, Pastel2, Pastel1

Use with semantic meaning. For example:

- [ ] `RdYlGn` - For quality metrics (higher is better)
- [ ] `Blues`, `Greens`, `Purples`, `Oranges` - For neutral counts
- [ ] `RdYlBu` - For metrics where middle is best (e.g., GC%)
- [ ] `Reds`, `OrRd` - For error/warning counts

## Phase 5: Registration

### search_patterns.yaml

```yaml
toolname/subtool:
  contents_re: "^expected_header_pattern"
  num_lines: 1
```

- [ ] Use `fn` file filename matching where possible
- [ ] Use `contents` for exact match or `contents_re` for regex
- [ ] Set appropriate `num_lines` to check
- [ ] Place alphabetically in file

### pyproject.toml

```toml
toolname = "multiqc.modules.toolname:MultiqcModule"
```

- [ ] Add entry point alphabetically in `[project.entry-points."multiqc.modules.v1"]`

### **init**.py

```python
from .toolname import MultiqcModule

__all__ = ["MultiqcModule"]
```

## Phase 6: Testing

### Running MultiQC

- [ ] Run MultiQC on the test data provided, from the MultiQC/test-data repository
- [ ] Use the `--strict` flag in the MultiQC command to find internal linting errors

### Unit Tests

- [ ] Create test file with sample data as string constants
- [ ] Test parser with full output format
- [ ] Test parser with minimal output format
- [ ] Test edge cases:
  - [ ] Empty file (header only)
  - [ ] Invalid format
  - [ ] Windows paths
  - [ ] Stdin input
  - [ ] Fallback sample name

### Integration Tests

- [ ] Ensure test data exists in `test-data` repository
- [ ] Run: `pytest tests/test_modules_run.py -k "toolname" -v`
- [ ] Verify both `test_all_modules` and `test_ignore_samples` pass

## Phase 7: Quality Checks

### Pre-commit

```bash
pre-commit run --files multiqc/modules/toolname/*
```

### Linting

```bash
ruff check multiqc/modules/toolname/
```

### Custom Checks

```bash
python .github/workflows/code_checks.py
```

Verifies:

- [ ] `add_data_source` is called
- [ ] `write_data_file` is called
- [ ] `doi=` is present in module init
- [ ] `add_software_version` is called

### Type Checking (if configured)

```bash
mypy multiqc/modules/toolname/
```

## Phase 8: Documentation

### Module Docstring

- [ ] Add comprehensive docstring to `MultiqcModule` class
- [ ] List all supported subcommands
- [ ] Document any configuration options
- [ ] Include example command to generate compatible output

### No Separate Markdown Files

- Do NOT create separate `.md` documentation files
- All documentation goes in the module class docstring

## Phase 9: Commit and PR

### Commit Message Format

```
Add [ToolName] module for [description]

Implements a new MultiQC module to parse output from [tool command].
[Tool description and purpose].

Features:
- [Feature 1]
- [Feature 2]

Closes #XXXX
```

### PR Description

- [ ] Reference original issue
- [ ] Describe metrics captured
- [ ] Note configuration options if any
- [ ] Reference test data location
- [ ] Include sample report screenshot if possible
