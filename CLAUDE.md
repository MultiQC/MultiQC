# CLAUDE.md

# IMPORTANT

Never push to main branch. When asked to add a change or create a pull request, always check if we are on main first. If we are, create a new branch and push to it.

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Testing

```bash
# Run all tests with coverage
pytest -vv --cov=multiqc --cov-report=xml

# Run tests in parallel (faster)
pytest -vv -n 4 --cov=multiqc --cov-report=xml

# Run single test file
pytest tests/test_modules_run.py -v

# Run specific test
pytest tests/test_modules_run.py::test_module_run -v
```

### Linting and Type Checking

```bash
# Run pre-commit hooks (includes ruff, mypy, etc.)
pre-commit run --all-files

# Run ruff linting
ruff check multiqc/

# Run mypy type checking
mypy multiqc
mypy tests

# Run custom code style checks
python .github/workflows/code_checks.py
```

### Installation and Development Setup

```bash
# Install MultiQC in development mode with all dependencies
pip install -e '.[dev]'

# Install just the package
pip install -e .
```

### Running MultiQC

```bash
# Basic usage
multiqc .

# Run with specific modules
multiqc . --module fastqc --module samtools

# Run with custom config
multiqc . --config custom_config.yaml
```

## High-Level Architecture

### Core Components

**Main Entry Point**: `multiqc/__main__.py` and `multiqc/multiqc.py`

- The `run()` function orchestrates the entire MultiQC execution
- Handles command-line parsing, config loading, and module execution

**Configuration System**: `multiqc/config.py`

- Global configuration loaded from `config_defaults.yaml`
- Supports user config files in multiple locations
- Environment variable support via `MULTIQC_*` variables
- Uses `pyproject.toml` for package configuration

**Module System**: `multiqc/modules/`

- Each bioinformatics tool has its own module directory
- Modules inherit from `BaseMultiqcModule` in `multiqc/base_module.py`
- Module discovery via entry points in `pyproject.toml`
- Pattern matching for log files defined in `search_patterns.yaml`

**Report Generation**: `multiqc/report.py`

- Collects data from all modules
- Generates HTML reports using Jinja2 templates
- Supports multiple output formats (HTML, data files, plots)

### Key Data Flow

1. **File Discovery**: `multiqc/core/file_search.py` scans directories for recognized log files
2. **Module Execution**: `multiqc/core/exec_modules.py` runs modules against found files
3. **Data Collection**: Modules parse log files and extract metrics
4. **Report Assembly**: `multiqc/core/write_results.py` combines all module outputs
5. **Output Generation**: Templates in `multiqc/templates/` render the final report

### Module Structure

All modules follow a consistent pattern:

- Inherit from `BaseMultiqcModule`
- Use `find_log_files()` to locate relevant files
- Parse log files and extract metrics
- Add data to general stats with `general_stats_addcols()`
- Create plots using classes in `multiqc/plots/`
- Add sections to the report with `add_section()`

### Plot System

**Plot Types**: `multiqc/plots/`

- Bar graphs, line graphs, scatter plots, heatmaps, tables
- Interactive plots using Plotly
- Flat plot export for publications

**Templates**: `multiqc/templates/`

- `default/`: Main HTML template with interactive features
- `simple/`: Simplified template
- `sections/`: Section-based template
- Template discovery via entry points

### AI Integration

**AI Features**: `multiqc/core/ai.py`

- Optional AI-powered report summaries
- Supports multiple providers (OpenAI, Anthropic, AWS Bedrock)
- Configurable via `ai_*` config options
- Sample name anonymization for privacy

## Important Files

- `multiqc/config_defaults.yaml`: Default configuration values
- `multiqc/search_patterns.yaml`: File pattern matching for all modules
- `pyproject.toml`: Package configuration and module entry points
- `multiqc/base_module.py`: Base class that all modules inherit from
- `multiqc/plots/`: Plot generation classes
- `multiqc/core/`: Core functionality (file search, module execution, etc.)

## Development Notes

- All modules are auto-discovered via entry points in `pyproject.toml`
- File patterns in `search_patterns.yaml` determine which files each module processes
- The codebase uses type hints extensively with mypy checking
- Tests are located in `tests/` and use pytest framework
- Code formatting is enforced with ruff
- The project supports Python 3.9+
- Use f-strings and other MODERN Python 3 syntax. Do not use `__future__` imports or `OrderedDict`'s.
- Use double quotes for strings.
- Do not add shebang lines to Python files unless they are placed in the `scripts/` folder.
- When writing modules, you must follow the following rules:
  - Raise `ModuleNoSamplesFound` when no samples are found. DO NOT RAISE `UserWarning`!
  - Call `self.add_software_version()`, even if version is not found, as it's required by linting.
  - Call `self.write_data_file` in the very end of the module, after all sections are added. IT IS IMPORTANT TO CALL IT IN THE END!
  - Add entry point into `pyproject.toml`. Ignore `setup.py`.
  - Do not add separate markdown files or module-level docstrings. Instead add a docstring to the module class.
  - Module's `info` MUST start with a capital letter.

DO NOT CREATE Pipfile - use pyproject.toml instead.
