"""
Generate docs/modules/<module.md> files for each module in MultiQC from the MultiqcModule class docstrings.

Usage:

    python scripts/generate_module_docs.py
    python scripts/generate_module_docs.py --output-dir docs/modules
    python scripts/generate_module_docs.py --module samtools
"""

import argparse
from textwrap import dedent

from multiqc import config, report, BaseMultiqcModule

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--module", help="Generate docs for a specific module")
args = parser.parse_args()


for mod_id, entry_point in config.avail_modules.items():
    if args.module and args.module != mod_id:
        continue
    if mod_id == "custom_content":
        continue

    output_path = config.REPO_DIR / "docs/modules" / f"{mod_id}.md"
    if output_path.exists():
        print(
            f"{output_path} already exists, so doing reverse: reading the md file contents, stripping the header, "
            f"and adding the rest to the module docstring, then writing the module code back to the python file to disk"
        )
        with output_path.open("r") as fh:
            text = fh.read()
        docstring = text.split("---\n", maxsplit=2)[2].strip()
        py_path = config.REPO_DIR / "multiqc" / "modules" / mod_id / f"{mod_id}.py"
        with py_path.open() as fh:
            module_code_lines = fh.readlines()
        with py_path.open("w") as fh:
            for line in module_code_lines:
                if "# Initialise the parent object" in line or "# Initialise the logger" in line:
                    continue

                fh.write(line)
                if line == "class MultiqcModule(BaseMultiqcModule):\n":
                    fh.write('    """\n')
                    for doc_line in docstring.split("\n"):
                        fh.write(f"    {doc_line}\n")
                    fh.write('    """\n')
                    fh.write("\n")

        print(f"Updated {py_path}")
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
url: {module.href}
description: {module.info}
---
{dedent(docstring)}
"""

    with output_path.open("w") as fh:
        fh.write(text)

    print(f"Generated {output_path}")
