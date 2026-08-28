# Module Implementation Checklist

Workflow guide for implementing a new MultiQC module. Hard requirements
(things that will fail lint or CI) are marked **must**; everything else is
guidance — adapt to the tool.

## Contents

- Research the tool
- Module skeleton (single-tool vs multi-subtool)
- Parser must-calls
- Error handling
- General stats columns
- Detailed sections and plots
- Section alerts
- Color scale guidelines
- Registration (search_patterns.yaml, pyproject.toml)
- Testing
- Quality checks
- Docstring and PR

## Research the tool

Read the tool's docs and any example output files. Note:

- Output filename / format (TSV, JSON, key-value, etc.)
- Whether the tool has subcommands with different output formats
- Version-specific variations and flag-dependent output (`--all`, `--verbose`)

Check `MultiQC/test-data/data/modules/` for existing test data; if absent,
ask for it or pull it from the request issue.

## Module skeleton

For a tool with one output format → single-tool module.
For a tool with multiple subcommands that emit different formats →
multi-subtool module (one parser file per subcommand).

See [module-structure.md](module-structure.md) for directory layout and
templates.

## Parser must-calls

Inside the parser loop, for each sample:

- **must** `module.add_data_source(f, s_name=s_name, section="...")`
- **must** `module.add_software_version(version_or_None, s_name)` — even
  when no version is available; lint requires it
- Use `module.clean_s_name(s_name, f)` only if `s_name` was extracted from
  file contents (not from `f["s_name"]`)
- If the tool has a standard extension that should be stripped, add it to
  `fn_clean_exts` in `config_defaults.yaml`

After the loop:

- **must** `data = module.ignore_samples(data)`
- **must** raise `ModuleNoSamplesFound` when `len(data) == 0` — never
  `UserWarning`
- **must** `module.write_data_file(data, "multiqc_toolname_subtool")` at
  the END, after all sections are added

For multi-subtool modules: the per-subtool `parse_*` function returns
`len(data)`. The orchestrator raises `ModuleNoSamplesFound` if every
subtool returned 0.

## Error handling

Don't silently default known fields. For output from trusted tools where a
value should always be present, access dict keys directly
(`parsed["total_reads"]`) instead of `parsed.get("total_reads", 0)`. A
missing key indicates a real format change — surfacing the `KeyError`
beats producing a report full of fake zeros. Reserve `.get(default)` for
genuinely optional fields. Catching a parse error to raise a friendlier
message with the file path is fine; the badness is silently producing fake
data.

For files that may not actually be from this tool (loose search pattern),
catch the parse error, `log.debug()` it, and `continue`.

## General stats columns

Keep the number of columns small (the general stats table is shared with
every other module). Per column, set:

- `title` — short display title
- `description` — tooltip text; use `config.read_count_desc` /
  `config.base_count_desc` etc. for shared formatting
- `scale` — see Color scale guidelines below
- `shared_key` — for related columns across modules (e.g. `"read_count"`)
- `hidden` — `True` for less-important columns
- `min` / `max` — for percentages and bounded metrics
- `suffix` — units (`"%"`, `" bp"`)

Only set keys that actually differ from defaults. Run through
`module.get_general_stats_headers()` for config integration, then
`module.general_stats_addcols(data, headers, namespace="toolname")`.

## Detailed sections and plots

`module.add_section()` arguments:

- `name` — section title, **human-readable** (`"Adapter Content"`, not
  `"adapter_content"`)
- `anchor` — unique section anchor
- `description` — section description
- `plot` — plot object (table, bargraph, linegraph, heatmap, violin, box,
  scatter)
- `alerts` — see Section alerts below

All user-facing strings (plot titles, axis labels `xlab`/`ylab`, table
column titles, section names) **must be human-readable**. Never pass raw
parsed keys like `total_counts` or `pct_dup` — convert to title case with
spaces.

Don't pre-filter samples in the main data dict at parse time. Keep every
sample in `self.toolname_data`; filter at plot-render time. This keeps
`write_data_file` complete and lets different sections make different
display choices.

## Section alerts

For warnings, notes, or "these samples were hidden" messages, use the
`alerts=` parameter on `add_section()` with a `SectionAlert` (from
`multiqc.types`):

- Pass affected sample names as `affected_samples=[...]` — they render in
  an expandable list automatically. **Don't** wrap names in `<code>`
  manually
- `level` must be a Bootstrap variant: `"primary"`, `"secondary"`,
  `"success"`, `"danger"`, `"warning"`, `"info"`, `"light"`, `"dark"`
- **Don't** append raw `<div class="alert ...">` HTML to `description` —
  that's the old pattern

If every sample is zero/empty for a plot, **keep the section visible**:
pass `plot=None` and add a `SectionAlert` so the user sees the analysis
ran and which samples were affected. Never drop the whole section.

See [code-patterns.md](code-patterns.md) for code examples.

## Color scale guidelines

Use ColorBrewer scales with semantic meaning:

- `RdYlGn` — quality metrics (higher is better)
- `RdYlBu` — middle-is-best (e.g. GC%)
- `Blues`, `Greens`, `Purples`, `Oranges` — neutral counts
- `Reds`, `OrRd` — error/warning counts

Available scales:

- Sequential: `OrRd`, `PuBu`, `BuPu`, `Oranges`, `BuGn`, `YlOrBr`, `YlGn`,
  `Reds`, `RdPu`, `Greens`, `YlGnBu`, `Purples`, `GnBu`, `Greys`,
  `YlOrRd`, `PuRd`, `Blues`, `PuBuGn`
- Diverging: `Spectral`, `RdYlGn`, `RdBu`, `PiYG`, `PRGn`, `RdYlBu`,
  `BrBG`, `RdGy`, `PuOr`
- Qualitative: `Set2`, `Accent`, `Set1`, `Set3`, `Dark2`, `Paired`,
  `Pastel2`, `Pastel1`

## Registration

`search_patterns.yaml` (alphabetical):

```yaml
toolname/subtool:
  contents_re: "^expected_header_pattern"
  num_lines: 1
```

Prefer `fn` (filename glob) when the tool produces a standard filename.
Fall back to `contents` for an exact-string match, or `contents_re` for a
regex. Set `num_lines` to cap content scanning; if the tool's stdout is
captured to a file, add 3 to the expected line count to allow for
prepended system lines.

**Audit upstream source** before finalising the pattern. Skim the tool's
code/docs for output variations across versions, optional flags, locale
settings. Matching one header line is fragile if the header changes
between versions. List version + flag combinations the pattern must
cover, then pick the smallest pattern that catches all without false
positives.

`pyproject.toml` — add alphabetically under
`[project.entry-points."multiqc.modules.v1"]`:

```toml
toolname = "multiqc.modules.toolname:MultiqcModule"
```

`__init__.py`:

```python
from .toolname import MultiqcModule

__all__ = ["MultiqcModule"]
```

## Testing

Run MultiQC against test data with `--strict` to surface internal lint
errors:

```bash
multiqc path/to/test-data/data/modules/toolname -m toolname --strict
```

Unit tests for each parser function. Cover:

- Full output format (all optional columns)
- Minimal output format
- Empty file (header only)
- Invalid format
- Windows paths
- Stdin input (fallback sample name)

Integration tests:

```bash
pytest tests/test_modules_run.py -k "toolname" -v
```

Both `test_all_modules` and `test_ignore_samples` must pass.

## Quality checks

```bash
prek run --files multiqc/modules/toolname/*
ruff check multiqc/modules/toolname/
python .github/workflows/code_checks.py
mypy multiqc/modules/toolname/
```

`code_checks.py` verifies `add_data_source`, `write_data_file`, `doi=`,
and `add_software_version` are present.

## Docstring and PR

The `MultiqcModule` class docstring is the module's user-facing
documentation. Include:

- What the tool does
- Supported subcommands (if multi-subtool)
- Any configuration options the module respects
- Example command(s) that generate compatible output

Do **not** add separate `.md` files for the module.

PR body: brief description first, then a `<details>` block with the full
write-up. Reference the original issue (`Closes #XXXX`), describe captured
metrics, mention any new config options, and link the test data location.
A sample report screenshot helps reviewers.
