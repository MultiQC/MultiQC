"""Shared schema + config_defaults.yaml loader used by both generator scripts."""

import sys
from pathlib import Path
from typing import AbstractSet, Any, Dict, List, Set, Tuple

import yaml

# Add repo root to path so scripts can import multiqc when run directly.
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import MultiQCConfig

DEFAULTS_PATH = Path(__file__).parent.parent / "multiqc" / "config_defaults.yaml"


def load_schema_and_defaults() -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
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


def load_sections(properties: Dict[str, Any], skip: AbstractSet[str] = frozenset()) -> Dict[str, List[str]]:
    """Group property names by their ``section`` tag, preserving source order.

    Sections appear in the order they first occur in ``properties`` (which
    matches the source order of ``with section(...)`` blocks in
    ``MultiQCConfig``). Within each section, field order also matches source
    order. Fails loudly if any property is missing a section tag and is not in
    ``skip``, so a new field can't be silently dropped from the docs/wizard.
    """
    sections: Dict[str, List[str]] = {}
    untagged: List[str] = []
    for prop_name, prop in properties.items():
        if prop_name in skip:
            continue
        section = prop.get("section")
        if section is None:
            untagged.append(prop_name)
            continue
        sections.setdefault(section, []).append(prop_name)
    if untagged:
        raise RuntimeError(
            f"{len(untagged)} schema property/properties have no 'section' tag: {sorted(untagged)}.\n"
            f'Wrap each Field with cfg(..., section="...") in multiqc/utils/config_schema.py, '
            f"or add to the loader caller's skip set."
        )
    return sections


def load_uncommon(properties: Dict[str, Any]) -> Set[str]:
    """Set of property names flagged ``advanced=True`` in their cfg() call."""
    return {name for name, prop in properties.items() if prop.get("advanced")}
