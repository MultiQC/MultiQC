"""
MultiQC configuration file validator script.
"""

import logging
import sys
import traceback
from pathlib import Path

import yaml
from pydantic import ValidationError
from rich.console import Console
from rich.panel import Panel
from rich.highlighter import ReprHighlighter

from multiqc.utils.config_schema import MultiQCConfig

logger = logging.getLogger(__name__)
console = Console()
highlighter = ReprHighlighter()


def check_config_file(config_file: str) -> None:
    """
    Check a MultiQC configuration file for errors.

    Args:
        config_file: Path to the configuration file to check
    """

    path = Path(config_file)
    if not path.exists():
        logger.error(f"Config file does not exist: {config_file}")
        sys.exit(1)

    if not path.is_file():
        logger.error(f"Config path is not a file: {config_file}")
        sys.exit(1)

    # Try to parse the file
    try:
        with open(config_file) as f:
            yaml_conf = yaml.safe_load(f)
            if yaml_conf is None:
                yaml_conf = {}  # Empty file
    except Exception as e:
        console.print(f"[bold red]YAML parsing error in config file:[/bold red] {str(e)}")
        console.print(traceback.format_exc())
        sys.exit(1)

    # Print parsed config for user to confirm
    console.print(f"[bold green]Successfully parsed YAML in {config_file}[/bold green]")

    # Validate against the Pydantic model
    has_validation_errors = False
    validation_errors = []

    try:
        # Validate using the MultiQCConfig model
        MultiQCConfig(**yaml_conf)
        console.print("[bold green]Config validation passed[/bold green]")

    except ValidationError as e:
        has_validation_errors = True
        validation_errors = e.errors()

        console.print("[bold red]Config validation failed:[/bold red]")
        for error in validation_errors:
            loc = ".".join(str(line) for line in error["loc"])
            console.print(f"[red]Error in [bold]{loc}[/bold]: {error['msg']}[/red]")

    if has_validation_errors:
        console.print("\n[bold red]Validation failed with errors.[/bold red]")
        sys.exit(1)
    else:
        console.print("\n[bold green]Config file passed all checks ✓[/bold green]")
        sys.exit(0)
