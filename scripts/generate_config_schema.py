#!/usr/bin/env python
"""
Script to generate JSON Schema for MultiQC config files.
"""

import json
from pathlib import Path

from multiqc.utils.config_schema import config_to_schema


def main():
    """Generate JSON schema file from Pydantic model"""
    schema = config_to_schema()

    schema_file = Path(__file__).parent.parent / "multiqc" / "utils" / "config_schema.json"
    with schema_file.open("w") as f:
        json.dump(schema, f, indent=4)

    print(f"Generated schema file: {schema_file}")


if __name__ == "__main__":
    main()
