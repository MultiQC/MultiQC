# MultiQC repository context

MultiQC is a Python tool that aggregates bioinformatics QC reports from many
upstream tools (FastQC, samtools, STAR, Salmon, etc.) into one HTML report.

## Categories that matter for triage

- **bug** — incorrect output, crashes, regressions
- **feature-request** — new functionality in core or templates
- **module-request** — support for a new bioinformatics tool. These usually
  arrive with the label `module: new` and are handled by a specialized agent.
- **docs** — documentation gaps
- **question** — usage support (often "how do I parse X")
- **other** — meta, governance, release coordination

## Areas of the codebase

- `multiqc/modules/<tool>/` — one directory per supported bioinformatics tool
- `multiqc/plots/` — plot generation (bar, line, scatter, heatmap, table)
- `multiqc/core/` — file search, module execution, report assembly
- `multiqc/templates/` — HTML report templates (default, simple, sections)
- `multiqc/config_defaults.yaml` — default configuration
- `multiqc/search_patterns.yaml` — file pattern matching per module

## Conventions reviewers should check

- Modules inherit from `BaseMultiqcModule`
- Raise `ModuleNoSamplesFound` when no samples found (never `UserWarning`)
- Call `self.add_software_version()` even when version is absent
- Call `self.write_data_file` at the end of the module
- Module `info` strings start with a capital letter
- Use f-strings, double quotes, modern Python 3 syntax (no `__future__`, no `OrderedDict`)
- Tests live in `tests/` (pytest)
