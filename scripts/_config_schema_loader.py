"""Shared schema + config_defaults.yaml loader used by both generator scripts."""

import sys
from pathlib import Path

import yaml

# Add repo root to path so scripts can import multiqc when run directly.
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import MultiQCConfig

DEFAULTS_PATH = Path(__file__).parent.parent / "multiqc" / "config_defaults.yaml"


def load_schema_and_defaults():
    """Load the JSON schema and config_defaults.yaml.

    Returns ``(properties, defaults, schema)``: the top-level ``properties`` dict,
    the parsed defaults YAML, and the full JSON schema (kept for callers that
    need to walk ``$defs``). Exits with a clear error if either source fails.
    """
    schema = MultiQCConfig.model_json_schema()
    try:
        with open(DEFAULTS_PATH, "r") as f:
            defaults = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {DEFAULTS_PATH}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)
    return schema.get("properties", {}), defaults, schema
