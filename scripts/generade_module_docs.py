"""
Generate docs/modules/<module.md> files for each module in MultiQC from the MultiqcModule class docstrings.

Usage:

    python scripts/generate_module_docs.py
    python scripts/generate_module_docs.py samtools
"""

import argparse
from textwrap import dedent

from multiqc import config, report, BaseMultiqcModule

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("module", help="Generate for a specific module", nargs="?")
args = parser.parse_args()


for mod_id, entry_point in config.avail_modules.items():
    if args.module and args.module != mod_id:
        continue
    if mod_id == "custom_content":
        continue

    data_dir = config.REPO_DIR / "test-data" / "data"
    mod_dir = data_dir / "modules" / mod_id
    assert mod_dir.exists() and mod_dir.is_dir()
    report.analysis_files = [mod_dir]
    report.search_files([mod_id])

    module_cls = entry_point.load()
    docstring = module_cls.__doc__
    if not docstring:
        print(f"Skipping {mod_id}: no docstring")
        continue

    module: BaseMultiqcModule = module_cls()

    text = f"""\
---
name: {module.name}
url: {module.href[0] if len(module.href) == 1 else module.href}
description: {module.info}
---
{dedent(docstring)}
"""

    output_path = config.REPO_DIR / "docs/modules" / f"{mod_id}.md"
    with output_path.open("w") as fh:
        fh.write(text)

    print(f"Generated {output_path}")
