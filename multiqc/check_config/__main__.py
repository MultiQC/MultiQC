#!/usr/bin/env python

"""
CLI helper tool to validate a MultiQC config file using the JSON schema.

$ multiqc-check multiqc_config.yml
"""

import yaml
import rich_click as click
from pathlib import Path
from multiqc.utils.config_schema import MultiQCConfig
from pydantic_core._pydantic_core import ValidationError
from rich.console import Console
from rich.panel import Panel
from rich.highlighter import ReprHighlighter
import sys

console = Console()
highlighter = ReprHighlighter()


@click.command()
@click.argument("config_path")
def run_multiqc_check(config_path):
    config_path = Path(config_path)
    with config_path.open() as f:
        config = yaml.safe_load(f)
    try:
        MultiQCConfig(**config)
    except ValidationError as e:
        console.print(f"\nError: MultiQC configuration file is invalid: [cyan]{config_path}\n", style="red")
        console.print(Panel(highlighter(str(e)), title="Validation Error", title_align="left", border_style="red"))
        console.print("")
        sys.exit(1)
    else:
        console.print(f"\nSuccess: MultiQC configuration file is valid: [cyan]{config_path}\n", style="green")


# Script is run directly
# NB: Usually runs with pyproject.toml console_scripts: multiqc-check = "multiqc.check_config.__main__:run_multiqc_check"
if __name__ == "__main__":
    run_multiqc_check()
