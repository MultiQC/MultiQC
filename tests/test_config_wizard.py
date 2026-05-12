"""Tests for scripts/generate_config_wizard.py."""

import importlib.util
import sys
from pathlib import Path

import yaml

from multiqc.utils.config_schema import MultiQCConfig

REPO_ROOT = Path(__file__).parent.parent
WIZARD_SCRIPT = REPO_ROOT / "scripts" / "generate_config_wizard.py"
CONFIG_DEFAULTS = REPO_ROOT / "multiqc" / "config_defaults.yaml"


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
