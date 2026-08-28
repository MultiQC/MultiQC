#!/usr/bin/env python
"""
Script to generate JSON Schema for MultiQC config files.
"""

import json
import sys
from pathlib import Path

# Prefer the repo's source tree over any installed `multiqc` in site-packages
# so this script always reflects local edits, with or without `pip install -e`.
sys.path.insert(0, str(Path(__file__).parent.parent))

from multiqc.utils.config_schema import config_to_schema  # noqa: E402


def main():
    """Generate JSON schema file from Pydantic model"""
    schema = config_to_schema()

    schema_file = Path(__file__).parent.parent / "multiqc" / "utils" / "config_schema.json"
    with schema_file.open("w") as f:
        # 2-space indent matches Prettier's default so the file survives the
        # commit hook without further edits.
        json.dump(schema, f, indent=2)
        f.write("\n")

    print(f"Generated schema file: {schema_file}")


if __name__ == "__main__":
    main()
