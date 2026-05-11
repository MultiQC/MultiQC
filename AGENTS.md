# AGENTS.md

## Cursor Cloud specific instructions

MultiQC is a standalone Python CLI tool with no external services required. See `CLAUDE.md` for development commands (install, test, lint, run).

### Key setup notes

- Python dev install: `pip install -e '.[dev]'`
- The `multiqc` binary installs to `~/.local/bin/` — ensure `$HOME/.local/bin` is on `PATH`.
- Tests require the test-data repo: `git clone --depth 1 https://github.com/MultiQC/test-data` in the workspace root. This is a separate repo and is gitignored.
- Frontend assets (JS/CSS) are built with Vite in `multiqc/templates/default/`. Run `npm install && npm run build` there if modifying frontend source files in `src/`.
- The `compiled/` directory under `multiqc/templates/default/` is auto-generated — never edit it directly.

### Running tests

```bash
pytest tests/ -v                    # all tests (~500, ~90s)
pytest tests/ -v -n 4              # parallel (faster)
pytest tests/test_cmdline.py -v    # single file
```

### Linting

```bash
ruff check multiqc/     # lint
mypy multiqc            # type check
```

### Running the application

```bash
multiqc test-data/data/modules/fastqc -o /tmp/multiqc_output --force
```

This generates a self-contained HTML report at the output path.
