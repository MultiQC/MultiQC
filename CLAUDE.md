# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## IMPORTANT

Never push to main. When asked to add a change or open a PR, check the
current branch first; if it's `main`, create a new branch and push to that.

Do NOT create a `Pipfile` — this project uses `pyproject.toml`.

## Style and conventions

- Python 3.9+. Use f-strings and other modern Python 3 syntax. Do not use
  `__future__` imports or `OrderedDict`.
- Double quotes for strings.
- No shebang lines on Python files unless under `scripts/`.
- Type hints are used extensively; mypy enforces.
- Code formatting is enforced with ruff.
- Tests live in `tests/` and use pytest.
- Helpers for distinct parts of the report — one method per section, one
  for parsing, one for general stats — are encouraged. They make
  `__init__` readable as a high-level outline. What to avoid is trivial
  wrappers: a helper that wraps one or two lines, just renames a
  one-liner, or is called once with no logical separation. Ask whether
  the function name is more meaningful than the code it hides; if not,
  inline it.
- Crash loudly on unexpected data. When parsing output from known
  bioinformatics tools, don't silently default known fields to empty
  dicts or zero values — that hides real format breakage behind a
  fake-looking report. Access documented fields directly
  (`parsed["total_reads"]`), not via `.get(key, 0)`. Catching a parse
  error to raise a friendlier message with the file path is fine;
  silently producing fake data is not. Reserve `.get(default)` for
  genuinely optional fields.
- Never use em-dashes (—) in any user-facing text: module descriptions,
  section titles, plot help text, docstrings, PR descriptions, commit
  messages. Use commas, semicolons, parentheses, or split into two
  sentences.
- Documentation examples in generic docs (e.g.
  `docs/markdown/development/modules.md`) should be tool-agnostic. Use
  placeholder names like `toolname` rather than referencing the specific
  module that motivated the example.

## Module rules (mandatory)

When writing modules, the following are mandatory:

- Raise `ModuleNoSamplesFound` when no samples are found. **Do not raise
  `UserWarning`**.
- Call `self.add_software_version()` even if version is unknown — it's
  required by linting.
- Call `self.write_data_file()` at the **very end** of the module, after
  all sections are added.
- Register the module via the entry point in `pyproject.toml` (ignore
  `setup.py`).
- Put module documentation in the module class docstring; do not add
  separate markdown files or module-level docstrings.
- The module's `info` field must start with a capital letter.

For full module guidance, see the `implementing-new-modules` skill in
`.claude/skills/`.

## Architecture

**Entry point**: `multiqc/__main__.py` and `multiqc/multiqc.py`. The
`run()` function orchestrates execution: CLI parsing, config loading,
module execution, report writing.

**Config**: `multiqc/config.py` loads defaults from
`multiqc/config_defaults.yaml` and supports user configs and `MULTIQC_*`
env vars.

**Modules**: `multiqc/modules/<toolname>/`. Each module inherits from
`BaseMultiqcModule` (`multiqc/base_module.py`). Modules are
auto-discovered via entry points in `pyproject.toml`. Search patterns
live in `multiqc/search_patterns.yaml`.

**Plots**: `multiqc/plots/` — bar graphs, line graphs, scatter,
heatmaps, tables, etc. Interactive plots use Plotly; flat plots are
also supported.

**Report**: `multiqc/report.py` collects module outputs and renders via
Jinja2 templates in `multiqc/templates/` (default, simple, sections,
gathered, original, geo, disco).

**Data flow**:

1. `multiqc/core/file_search.py` scans for recognised log files.
2. `multiqc/core/exec_modules.py` runs modules against found files.
3. Modules parse and emit data.
4. `multiqc/core/write_results.py` assembles the final report.

**AI features** (optional): `multiqc/core/ai.py`. Multi-provider
(OpenAI, Anthropic, AWS Bedrock), config via `ai_*` options. Supports
sample name anonymisation.

## Key files

- `multiqc/config_defaults.yaml` — default configuration values
- `multiqc/utils/config_schema.py` — Pydantic `MultiQCConfig` model: source
  of truth for option types, descriptions, examples, and Literal enums.
  Driven into `docs/markdown/config_schema.md`,
  `multiqc/utils/config_schema.json`, and the wizard HTML by scripts in
  `scripts/`. `tests/test_config_wizard.py` catches drift.
- `multiqc/search_patterns.yaml` — file-pattern matching for every module
- `multiqc/base_module.py` — `BaseMultiqcModule` parent class
- `pyproject.toml` — package config and module entry points
- `multiqc/plots/` — plot classes
- `multiqc/core/` — file search, module execution, write results

## Development commands

```bash
# Tests
pytest -vv --cov=multiqc --cov-report=xml          # full suite
pytest -vv -n 4 --cov=multiqc --cov-report=xml     # parallel
pytest tests/test_modules_run.py -v                # one file
pytest tests/test_modules_run.py::test_module_run -v   # one test

# Lint / type-check
prek run --all-files
ruff check multiqc/
mypy multiqc
mypy tests
python .github/workflows/code_checks.py            # custom checks

# Install for development
pip install -e '.[dev]'

# Run multiqc
multiqc .                                          # all modules
multiqc . --module fastqc --module samtools        # specific
multiqc . --config custom_config.yaml              # custom config
```
