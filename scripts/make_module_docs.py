"""
Generate docs/modules/<module.md> files for each module in MultiQC from the MultiqcModule class docstrings.

Usage:

    python scripts/make_module_docs.py
    python scripts/make_module_docs.py samtools
"""

import argparse
from pathlib import Path
from textwrap import dedent, indent

from multiqc import config, report, BaseMultiqcModule

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("module", help="Generate for a specific module", nargs="?")
args = parser.parse_args()


for mod_id, entry_point in config.avail_modules.items():
    if args.module and args.module != mod_id:
        continue
    if mod_id == "custom_content":
        continue

    mod_dir = Path("test-data") / "data" / "modules" / mod_id
    assert mod_dir.exists() and mod_dir.is_dir(), mod_dir
    report.analysis_files = [mod_dir]
    report.search_files([mod_id])

    module_cls = entry_point.load()
    docstring = module_cls.__doc__ or ""

    module: BaseMultiqcModule = module_cls()

    if module.extra:
        extra = "\n".join(line.strip() for line in module.extra.split("\n") if line.strip())
        extra = "\nextra_description: >\n" + indent(extra, "  ")
    else:
        extra = ""
    text = f"""\
---
name: {module.name}
urls: {module.href}
summary: >
  {module.info}{extra}
---

{dedent(docstring)}
""".strip()

    # Remove double empty lines
    while "\n\n\n" in text:
        text = text.replace("\n\n\n", "\n\n")

    output_path = Path("docs") / "modules" / f"{mod_id}.md"
    with output_path.open("w") as fh:
        fh.write(text)

    print(f"Generated {output_path}")
