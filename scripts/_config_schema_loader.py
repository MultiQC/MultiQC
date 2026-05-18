"""Shared schema + config_defaults.yaml loader used by both generator scripts."""

import sys
from pathlib import Path
from typing import AbstractSet, Any, Dict, List, Optional, Tuple

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


def load_sections_with_groups(
    properties: Dict[str, Any],
    skip: AbstractSet[str] = frozenset(),
) -> Dict[str, Dict[Optional[str], List[str]]]:
    """Group property names by ``section`` then by ``group``, preserving source order.

    Returns a dict keyed by section name, where each value is a dict keyed by
    group name (or ``None`` for ungrouped fields). Sections, groups within a
    section, and fields within a group all appear in source order. Group
    members must be contiguous in source order; non-contiguous groups produce
    duplicate headings in the generated docs and wizard. Fails loudly on any
    property missing a ``section`` tag that is not in ``skip``.
    """
    out: Dict[str, Dict[Optional[str], List[str]]] = {}
    untagged: List[str] = []
    for prop_name, prop in properties.items():
        if prop_name in skip:
            continue
        section = prop.get("section")
        if section is None:
            untagged.append(prop_name)
            continue
        group = prop.get("group")
        section_buckets = out.setdefault(section, {})
        section_buckets.setdefault(group, []).append(prop_name)
    if untagged:
        raise RuntimeError(
            f"{len(untagged)} schema property/properties have no 'section' tag: {sorted(untagged)}.\n"
            f'Wrap each Field with cfg(..., section="...") in multiqc/utils/config_schema.py, '
            f"or add to the loader caller's skip set."
        )
    return out
