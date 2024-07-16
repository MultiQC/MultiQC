"""
Usage:

    python scripts/doc_md_to_docstrings.py
    python scripts/doc_md_to_docstrings.py samtools
"""

import argparse

import yaml

from multiqc import config

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("module", help="Generate for a specific module", nargs="?")
args = parser.parse_args()


for mod_id, entry_point in config.avail_modules.items():
    if args.module and args.module != mod_id:
        continue
    if mod_id == "custom_content":
        continue

    output_path = config.REPO_DIR / "docs/modules" / f"{mod_id}.md"
    if not output_path.exists():
        print(f"md file {output_path} doesn't exist, skipping")
        continue

    with output_path.open("r") as fh:
        text = fh.read()

    descr = yaml.safe_load(text.split("---\n")[1].strip())["description"]
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
                fh.write(f"    {descr}\n")
                fh.write("\n")
                for doc_line in docstring.split("\n"):
                    if doc_line.strip():
                        fh.write(f"    {doc_line}\n")
                    else:
                        fh.write("\n")  # ruff would strip spaces from empty lines, so why not do this before him
                fh.write('    """\n')
                fh.write("\n")

    print(f"Updated {py_path}")
