"""Tests for scripts/generate_config_wizard.py."""

import importlib.util
import keyword
import re
import sys
from pathlib import Path

import yaml

from multiqc.utils.config_schema import MultiQCConfig

REPO_ROOT = Path(__file__).parent.parent
WIZARD_SCRIPT = REPO_ROOT / "scripts" / "generate_config_wizard.py"
CONFIG_DEFAULTS = REPO_ROOT / "multiqc" / "config_defaults.yaml"
CONFIG_PY = REPO_ROOT / "multiqc" / "config.py"

# Top-level annotated names in multiqc/config.py that are runtime / derived
# state, not user-settable YAML keys. The schema regression test below uses
# this to skip non-config attributes. Add new runtime fields here; add new
# user-facing fields to the schema itself.
RUNTIME_CONFIG_ATTRS = {
    "analysis_dir",
    "avail_modules",
    "avail_templates",
    "data_dir",
    "explicit_user_config_files",
    "filename",
    "kwargs",
    "loaded_user_files",
    "modules_dir",
    "nondefault_config",
    "output_dir",
    "output_fn",
    "plots_dir",
    "working_dir",
}


def _load_wizard_module():
    spec = importlib.util.spec_from_file_location("generate_config_wizard", WIZARD_SCRIPT)
    assert spec is not None and spec.loader is not None, f"Could not load {WIZARD_SCRIPT}"
    module = importlib.util.module_from_spec(spec)
    sys.modules["generate_config_wizard"] = module
    spec.loader.exec_module(module)
    return module


def test_every_field_has_a_section():
    """Every MultiQCConfig field must declare a section via cfg(..., section=...).

    Section tags drive both the config docs grouping and the wizard sidebar.
    A field with no section would be silently dropped from both, so the
    generator scripts refuse to build until one is set. Add the section to the
    Field via cfg() in `multiqc/utils/config_schema.py`, or, if it cannot
    reasonably be rendered, add it to `SKIP_PROPERTIES` in the generator script.
    """
    wizard = _load_wizard_module()
    properties = MultiQCConfig.model_json_schema()["properties"]
    untagged = sorted(
        name for name, prop in properties.items() if name not in wizard.SKIP_PROPERTIES and "section" not in prop
    )
    assert not untagged, (
        f"Config properties with no section tag: {untagged}. "
        f'Wrap each Field with cfg(..., section="...") in multiqc/utils/config_schema.py.'
    )


def test_wizard_skip_list_is_in_schema():
    """SKIP_PROPERTIES must reference real config fields, not stale names."""
    wizard = _load_wizard_module()
    schema_props = set(MultiQCConfig.model_json_schema()["properties"])
    stale = wizard.SKIP_PROPERTIES - schema_props
    assert not stale, f"SKIP_PROPERTIES references unknown fields: {sorted(stale)}"


def test_config_py_annotations_are_in_schema():
    """Every user-facing annotation in multiqc/config.py must have a matching
    schema field. Without this guard a real-world YAML config (eg. from a
    pipeline) flags valid options as "unknown" in the wizard's Validate view
    and Pydantic can't validate them at config-load time. Add runtime-only
    fields (derived paths, entry-point lookups, internal state) to
    ``RUNTIME_CONFIG_ATTRS`` in this file rather than to the schema.
    """
    attrs = set()
    with open(CONFIG_PY) as f:
        for line in f:
            m = re.match(r"^([a-zA-Z][a-zA-Z0-9_]*)\s*:\s", line)
            if m and not keyword.iskeyword(m.group(1)):
                attrs.add(m.group(1))
    schema_props = set(MultiQCConfig.model_json_schema()["properties"])
    user_facing = attrs - RUNTIME_CONFIG_ATTRS
    missing = sorted(user_facing - schema_props)
    assert not missing, (
        f"multiqc/config.py annotates these user-facing keys that are missing from MultiQCConfig: {missing}. "
        f"Add a Field for each in multiqc/utils/config_schema.py, "
        f"or — if it's runtime-only — list it in RUNTIME_CONFIG_ATTRS in this test."
    )


def test_runtime_config_attrs_still_exist():
    """Stale names in RUNTIME_CONFIG_ATTRS hide drift from config.py: if a
    runtime attribute is renamed or removed there, this test forces the
    allow-list to follow suit.
    """
    attrs = set()
    with open(CONFIG_PY) as f:
        for line in f:
            m = re.match(r"^([a-zA-Z][a-zA-Z0-9_]*)\s*:\s", line)
            if m and not keyword.iskeyword(m.group(1)):
                attrs.add(m.group(1))
    stale = sorted(RUNTIME_CONFIG_ATTRS - attrs)
    assert not stale, (
        f"RUNTIME_CONFIG_ATTRS references names no longer in multiqc/config.py: {stale}. "
        f"Remove them from the set in tests/test_config_wizard.py."
    )


def test_config_defaults_keys_are_in_schema():
    """Every key in config_defaults.yaml must have a matching field in MultiQCConfig.

    Without this guard, a developer can add a config option (with a default in the
    YAML and a typed attribute in multiqc/config.py) but forget to surface it in
    the schema, leaving it undocumented and absent from the wizard.
    """
    with open(CONFIG_DEFAULTS) as f:
        defaults = yaml.safe_load(f) or {}
    schema_props = set(MultiQCConfig.model_json_schema()["properties"])
    missing = sorted(set(defaults) - schema_props)
    assert not missing, (
        f"Config defaults present in config_defaults.yaml but missing from MultiQCConfig: {missing}. "
        f"Add a Field for each to multiqc/utils/config_schema.py."
    )
