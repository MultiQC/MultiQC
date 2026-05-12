---
name: implementing-new-modules
description: |
  Create a new MultiQC module from scratch. Parse a bioinformatics tool's output, register search patterns and entry points, add general stats columns, build plots and sections, write tests, open a PR. Use when implementing a `module: new` GitHub issue, when the user asks to add support for a new tool, or when adding a parser for a new tool output format.
---

# Implement New MultiQC Module

## Workflow

1. **Research the tool**: read its docs, get example output files (check
   `MultiQC/test-data/data/modules/` first), find similar existing
   modules for reference. Note version- and flag-dependent output
   variations.
2. **Pick architecture**: single-tool or multi-subtool (see below).
3. **Build the module**: parser, general stats, sections, plots.
   Full step-by-step in [implementation-checklist.md](implementation-checklist.md).
4. **Register**: add to `search_patterns.yaml` and `pyproject.toml`
   entry points.
5. **Test**: unit tests for parsers, integration test via
   `pytest tests/test_modules_run.py -k "toolname" -v`, and
   `multiqc … --strict` on test data.
6. **Quality gate**: `prek run`, `ruff check`, `python .github/workflows/code_checks.py`, `mypy`.
7. **PR**: brief summary, `<details>` for the full write-up, `Closes #XXXX`.

## Architecture decision

**Single-tool** (FastQC, Qualimap): one output format, one parser.

```
multiqc/modules/toolname/
├── __init__.py
└── toolname.py
```

**Multi-subtool** (samtools, seqkit, picard): distinct subcommands with
different output formats, or more subcommands likely to be added.

```
multiqc/modules/toolname/
├── __init__.py
├── toolname.py        # Orchestrator
├── subtool1.py        # parse_toolname_subtool1() function
├── subtool2.py        # parse_toolname_subtool2() function
└── tests/
```

Class skeletons and full file templates: [module-structure.md](module-structure.md).
Patterns for code _inside_ a module (parsing, plots, alerts, etc.):
[code-patterns.md](code-patterns.md).

## Common Pitfalls

1. **Forgetting `add_software_version()`** — required by linting, even if version is `None`.
2. **Calling `write_data_file()` too early** — must be at end, after all sections.
3. **Raising `UserWarning` instead of `ModuleNoSamplesFound`**.
4. **Not handling both tab- and space-separated output** when both are valid.
5. **Hardcoding values instead of using `f["s_name"]`** and other dynamic variables.
6. **Manually cleaning sample names instead of `self.clean_s_name()`**.
7. **Inappropriate colour scales** — e.g. `RdYlGn` for GC% (which is not "higher is better").
8. **Silently defaulting on known fields** — `parsed.get(key, 0)` for fields the tool always emits (or `try/except: return {}` over the whole parse) hides real format breakage behind a fake-looking report. Access documented keys directly; reserve `.get(default)` for genuinely optional fields. Catching a parse error to raise a friendlier message is fine; silently producing zeros is not.
9. **Trivial single-statement helpers** — a helper that wraps one or two lines, or just renames a one-liner, adds indirection without aiding readability. Per-section / per-parser helpers (`_add_adapter_section`, `_parse_log`) are fine and often clearer; one-liner wrappers (`_add_filtered_section` calling `add_section` with no real logic) are not.
10. **Using raw parsed dict keys in user-facing text** — `total_counts` and `pct_dup` belong in code, never in plot/column titles, axis labels, or section names. Convert to `"Total Counts"`, `"% Duplicates"`.
11. **Dropping the whole section when all samples are zero** — keep the section, pass `plot=None`, and add a `SectionAlert` via the `alerts=` parameter on `add_section()` listing affected samples. Don't append raw `<div class="alert ...">` HTML to `description` — use the `alerts` API.
12. **Pre-filtering samples at parse time** — keep every sample in the main data dict so `write_data_file` is complete; filter at plot-render time.
13. **Em-dashes (—) in any user-facing text** — descriptions, docstrings, alerts, PR text. AI tell. Use commas or split sentences.
