"""
MultiQC configuration file validator script.
"""

import logging
import sys
from pathlib import Path

import yaml
from pydantic import ValidationError
from rich.console import Console, group
from rich.highlighter import ReprHighlighter
from rich.panel import Panel
from rich.text import Text

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
        if hasattr(e, "problem") and hasattr(e, "problem_mark"):
            parsing_err = highlighter(f"{str(e.problem).strip()} {str(e.problem_mark).strip()}")
        else:
            parsing_err = Text("Could not parse YAML")
        console.print(
            Panel(
                parsing_err,
                title=f":x: Error - Could not parse configuration file YAML: [bold]{config_file}",
                title_align="left",
                border_style="red",
            )
        )
        sys.exit(1)

    # Print parsed config for user to confirm
    logger.debug(f"Successfully parsed YAML in {config_file}")

    try:
        # Validate using the MultiQCConfig model
        MultiQCConfig(**yaml_conf)
        console.print(f"[green]✅ Config validation passed: [bold]{config_file}")
        sys.exit(0)

    except ValidationError as e:

        @group()
        def get_error_msgs(e):
            validation_errors = e.errors()
            for error in validation_errors:
                loc = ".".join(str(line) for line in error["loc"])
                yield Text.from_markup(
                    f"[bold]Error in [blue]'{loc}'[/blue]:[/] {error['msg']} "
                    f"[dim](got: [green italic]{error['input'][:50]}[/green italic])[/]"
                )

        console.print(
            Panel(
                get_error_msgs(e),
                title=f":x: Error - MultiQC configuration file is invalid: [bold]{config_file}",
                title_align="left",
                border_style="red",
            )
        )
        sys.exit(1)
