# Implement New MultiQC Module

## Description

Create a new MultiQC module from a GitHub issue or feature request. This skill guides you through the complete process of implementing, testing, and submitting a new module for parsing bioinformatics tool output.

**Use this skill when:**

- Implementing a module from a `module: new` issue
- User asks to create a module for a specific tool
- Adding support for a new bioinformatics tool output format

## Prerequisites

Before starting implementation:

1. **Understand the tool output format** - Get sample output files or documentation
2. **Check test-data repository** - Look for existing test data in `MultiQC/test-data`
3. **Review similar modules** - Find modules for similar tools as reference

## Quick Start

1. **Research the tool** - Understand output format, check for subtools
2. **Create module structure** - Follow multi-subtool pattern if applicable
3. **Implement parser** - Parse output and extract metrics
4. **Add visualizations** - Tables, bar graphs, line plots as appropriate
5. **Register module** - Add entry point and search patterns
6. **Write tests** - Unit tests for parser, integration tests
7. **Run checks** - Linting, type checking, tests
8. **Submit PR** - With clear description and test data reference

## Module Architecture Decision

### Single-tool modules

For tools with one output format (e.g., FastQC, Qualimap):

```
multiqc/modules/toolname/
├── __init__.py
└── toolname.py        # MultiqcModule class with all logic
```

### Multi-subtool modules

For tools with multiple subcommands (e.g., samtools, seqkit, picard):

```
multiqc/modules/toolname/
├── __init__.py        # Exports MultiqcModule
├── toolname.py        # Orchestrator calling submodules
├── subtool1.py        # parse_toolname_subtool1() function
├── subtool2.py        # parse_toolname_subtool2() function
└── tests/
    ├── __init__.py
    └── test_subtool1.py
```

**Use multi-subtool pattern when:**

- Tool has distinct subcommands (e.g., `samtools stats`, `samtools flagstat`)
- Each subcommand has different output format
- Future subcommands are likely to be added

## Implementation Checklist

See [implementation-checklist.md](implementation-checklist.md) for detailed steps.

### Essential Steps

- [ ] Create module directory structure
- [ ] Implement parser function(s)
- [ ] Add to `search_patterns.yaml`
- [ ] Add entry point in `pyproject.toml`
- [ ] Call `self.add_data_source()`
- [ ] Call `self.add_software_version()`
- [ ] Call `self.write_data_file()` at end
- [ ] Raise `ModuleNoSamplesFound` when no data found
- [ ] Add general stats columns if needed
- [ ] Create detailed section(s) with table/plots
- [ ] Write unit tests

### Quality Checks

- [ ] `multiqc -f --strict path/to/test/data/`
- [ ] `pre-commit run`
- [ ] `ruff check multiqc/modules/yourmodule/`
- [ ] `python .github/workflows/code_checks.py`
- [ ] `pytest multiqc/modules/yourmodule/tests/ -v`
- [ ] `pytest tests/test_modules_run.py -k "yourmodule" -v`

## Code Patterns

See [code-patterns.md](code-patterns.md) for common patterns including:

- File discovery and parsing
- Sample name handling
- General stats headers
- Table and plot creation
- Multi-format parsing (tab/space-separated)

## File Registration

### search_patterns.yaml

If the tool has a standardised output filename, search by this as it's fastest:

```yaml
toolname/subtool:
  fn: "*_custom_extension.txt"
```

If this is not possible, search the file contents.
Use `contents` with a single string if possible, or the more expensive `contents_re`
if a regular expression is required.

Specify `num_lines` when searching by content to cap the amount of lines searched in the file.
If the tool saves stdout, add 3 to the number required, in case of the system prepending to the file.

```yaml
toolname/subtool:
  contents_re: "^header_pattern\\s+column1\\s+column2"
  num_lines: 1
```

### pyproject.toml

```toml
[project.entry-points."multiqc.modules.v1"]
toolname = "multiqc.modules.toolname:MultiqcModule"
```

## Testing Strategy

1. **Unit tests** - Test parser functions directly with sample data
2. **Integration tests** - Ensure module runs with test-data repository
3. **Edge cases** - Empty files, malformed input, stdin, Windows paths

## Common Pitfalls

1. **Forgetting `add_software_version()`** - Required by linting, even if version is None
2. **Calling `write_data_file()` too early** - Must be at end after all sections added
3. **Raising `UserWarning` instead of `ModuleNoSamplesFound`**
4. **Not handling both tab and space-separated output**
5. **Hardcoding values instead of using dynamic variables such as `f["s_name"]`**
6. **Manually cleaning sample names instead of using core functions like `self.clean_s_name()`**
7. **Using colour scales inappropriately, eg. RdYlGn scale for non-quality metrics** (GC% is not "higher is better")

## PR Submission

1. Reference the original issue: `Closes #XXXX`
2. Describe the tool and what metrics are captured
3. Note any configuration options added
4. Reference test data location in MultiQC/test-data
5. Include sample report screenshot if possible

Give the PR body a brief description, then add a detailed description of the work
done within a <details> tag.

## Files in This Skill

- `SKILL.md` (this file): Overview and workflow
- `implementation-checklist.md`: Detailed step-by-step checklist
- `code-patterns.md`: Common code patterns and examples
- `module-structure.md`: Directory and file structure templates

## Related Documentation

- [CLAUDE.md](../../../CLAUDE.md): Repository-specific instructions
- [Contributing Guide](../../../CONTRIBUTING.md): General contribution guidelines
- [Module Development Docs](https://docs.seqera.io/multiqc/development/modules): Official documentation
