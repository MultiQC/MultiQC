"""Shared schema + config_defaults.yaml loader used by both generator scripts."""

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

# Add repo root to path so scripts can import multiqc when run directly.
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import MultiQCConfig

MULTIQC_DIR = Path(__file__).parent.parent / "multiqc"
DEFAULTS_PATH = MULTIQC_DIR / "config_defaults.yaml"
SEARCH_PATTERNS_PATH = MULTIQC_DIR / "search_patterns.yaml"


def load_schema_and_defaults() -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Load the JSON schema and config_defaults.yaml.

    Returns ``(properties, defaults, schema)``: the top-level ``properties`` dict,
    the parsed defaults YAML, and the full JSON schema (kept for callers that
    need to walk ``$defs``). ``sp`` defaults are loaded from
    ``multiqc/search_patterns.yaml`` and folded into the defaults dict so
    callers see the same shape regardless of which file the value lives in
    at runtime. Exits with a clear error if any source fails.
    """
    schema = MultiQCConfig.model_json_schema()
    try:
        with open(DEFAULTS_PATH, "r") as f:
            defaults = yaml.safe_load(f) or {}
    except FileNotFoundError:
        print(f"Error: Could not find {DEFAULTS_PATH}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)
    try:
        with open(SEARCH_PATTERNS_PATH, "r") as f:
            defaults["sp"] = yaml.safe_load(f) or {}
    except FileNotFoundError:
        print(f"Error: Could not find {SEARCH_PATTERNS_PATH}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        sys.exit(1)
    return schema.get("properties", {}), defaults, schema


def load_sections_with_groups(
    properties: Dict[str, Any],
) -> Dict[str, Dict[Optional[str], List[str]]]:
    """Group property names by ``section`` then by ``group``, preserving source order.

    Returns a dict keyed by section name, where each value is a dict keyed by
    group name (or ``None`` for ungrouped fields). Sections, groups within a
    section, and fields within a group all appear in source order. Group
    members must be contiguous in source order; non-contiguous groups produce
    duplicate headings in the generated docs and wizard. Fails loudly on any
    property missing a ``section`` tag.
    """
    out: Dict[str, Dict[Optional[str], List[str]]] = {}
    untagged: List[str] = []
    for prop_name, prop in properties.items():
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
            f'Wrap each Field with cfg(..., section="...") in multiqc/utils/config_schema.py.'
        )
    return out
