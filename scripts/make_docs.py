"""
Generate docs/modules/<module.md> files for each module in MultiQC from the MultiqcModule class docstrings.

Usage:

<<<<<<< HEAD:scripts/make_module_docs.py
=======
    python scripts/make_module_docs.py
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py
    python scripts/make_module_docs.py modules
    python scripts/make_module_docs.py modules --module samtools
    python scripts/make_module_docs.py changelog
"""

import json
import os
from typing import Dict

import yaml

import argparse
from pathlib import Path
from textwrap import dedent, indent
import subprocess

from multiqc import config, report, BaseMultiqcModule

parser = argparse.ArgumentParser(description=__doc__)
<<<<<<< HEAD:scripts/make_module_docs.py
parser.add_argument("mode", help="Generate modules or changelog", nargs="?")
parser.add_argument("--module", help="Generate for a specific module", nargs="?")
args = parser.parse_args()
=======
parser.add_argument("--force", help="Force generation even if the files are up to date", action="store_true")
args = parser.parse_args()

DOCS_REPO = Path("seqeralabs-docs")
if not DOCS_REPO.exists():
    raise RuntimeError("Please clone https://github.com/seqeralabs/docs into seqeralabs-docs/")

TEST_DATA_DIR = Path("test-data")
if not TEST_DATA_DIR.exists():
    raise RuntimeError("Please clone https://github.com/MultiQC/test-data into test-data/")

OUTPUT_PATH = DOCS_REPO / "multiqc" / "docs"
os.makedirs(OUTPUT_PATH, exist_ok=True)

# Use rsync to sync the directories, except for README.md
rsync_command = [
    "rsync",
    "-av",  # archive mode, verbose
    "--delete",  # delete extraneous files from destination
    "docs/",  # source with trailing slash
    "--exclude",
    "README.md",
    f"{OUTPUT_PATH}/",  # destination with trailing slash
]

try:
    subprocess.run(rsync_command, check=True)
    print(f"Successfully synced docs/ to {OUTPUT_PATH}")
except subprocess.CalledProcessError as e:
    raise RuntimeError(f"Failed to sync directories: {e}")
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py

if args.mode and args.mode not in ("modules", "changelog"):
    parser.error("Invalid mode. Please use 'modules' or 'changelog'.")

if not args.mode or args.mode == "modules":
    test_data_dir = Path("test-data")
    if not test_data_dir.exists():
        raise RuntimeError("Please clone https://github.com/MultiQC/test-data into test-data/")

<<<<<<< HEAD:scripts/make_module_docs.py
    sp_by_mod: Dict[str, Dict] = dict()
    with (Path(config.MODULE_DIR) / "search_patterns.yaml").open() as f:
        for k, v in yaml.safe_load(f).items():
            mod_id = k.split("/")[0]
            sp_by_mod.setdefault(mod_id, {})[k] = v

    os.makedirs("docs/markdown/modules", exist_ok=True)
=======
# Table in the index page
modules_data = []
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py

    # Table in the index page
    modules_data = []

<<<<<<< HEAD:scripts/make_module_docs.py
    for mod_id, entry_point in config.avail_modules.items():
        if mod_id == "custom_content":
            continue
=======
    module_md_path = OUTPUT_PATH / "markdown" / "modules" / f"{mod_id}.md"
    if not args.force and module_md_path.exists():
        print(f"Skipping {module_md_path} because it already exists")
        continue

    mod_data_dir = TEST_DATA_DIR / "data" / "modules" / mod_id
    assert mod_data_dir.exists() and mod_data_dir.is_dir(), mod_data_dir
    report.analysis_files = [mod_data_dir]
    report.search_files([mod_id])
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py

        mod_data_dir = test_data_dir / "data" / "modules" / mod_id
        assert mod_data_dir.exists() and mod_data_dir.is_dir(), mod_data_dir
        report.analysis_files = [mod_data_dir]
        report.search_files([mod_id])

        module_cls = entry_point.load()
        docstring = module_cls.__doc__ or ""

<<<<<<< HEAD:scripts/make_module_docs.py
        module: BaseMultiqcModule = module_cls()
        modules_data.append(
            {"id": f"modules/{mod_id}", "data": {"name": f"{module.name}", "summary": f"{module.info}"}}
        )

        if args.module and args.module != mod_id:
            continue
=======
    if module.extra:
        extra = "\n".join(line.strip() for line in module.extra.split("\n") if line.strip())
        extra += "\n\n"
    else:
        extra = ""
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py

        if module.extra:
            extra = "\n".join(line.strip() for line in module.extra.split("\n") if line.strip())
            extra += "\n\n"
        else:
            extra = ""

        text = f"""\
---
title: {module.name}
displayed_sidebar: multiqcSidebar
description: >
<<<<<<< HEAD:scripts/make_module_docs.py
{module.info}
=======
{indent(module.info, "    ")}
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/{mod_id}/{mod_id}.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
{module.info}

{", ".join([f"[{href}]({href})" for href in module.href])}
:::

{extra}{dedent(docstring)}

### File search patterns

```yaml
{yaml.dump(sp_by_mod[mod_id]).strip()}
```
"""

        # Remove double empty lines
        while "\n\n\n" in text:
            text = text.replace("\n\n\n", "\n\n")

<<<<<<< HEAD:scripts/make_module_docs.py
        output_path = Path("docs") / "markdown" / "modules" / f"{mod_id}.md"
        with output_path.open("w") as fh:
            fh.write(text)

        print(f"Generated {output_path}")

    with (Path("docs") / "markdown" / "modules.mdx").open("w") as fh:
        fh.write(
            """\
=======
    module_md_path.parent.mkdir(parents=True, exist_ok=True)
    with module_md_path.open("w") as fh:
        fh.write(text)
    print(f"Generated {module_md_path}")

if not args.force and (mdx_path := OUTPUT_PATH / "markdown" / "modules.mdx").exists():
    print(f"Skipping {mdx_path} because it already exists")
else:
    with mdx_path.open("w") as fh:
        fh.write(
            f"""\
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py
---
title: Supported Tools
description: Tools supported by MultiQC
displayed_sidebar: multiqcSidebar
---
<<<<<<< HEAD:scripts/make_module_docs.py
            
=======
        
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py
MultiQC currently has modules to support {len(config.avail_modules)} bioinformatics tools, listed below.

Click the tool name to go to the MultiQC documentation for that tool.

:::tip[Missing something?]
If you would like another tool to to be supported, please [open an issue](https://github.com/MultiQC/MultiQC/issues/new?labels=module%3A+new&template=module-request.yml).
:::

import MultiqcModules from "@site/src/components/MultiqcModules";

<MultiqcModules
modules={{{str(json.dumps(modules_data))}}}
/>

"""
        )
<<<<<<< HEAD:scripts/make_module_docs.py


if not args.mode or args.mode == "changelog":
    """
    Script to sync the MultiQC changelog into the docs site.
    """

    # Read the changelog from the submodule
    changelog_path = Path("CHANGELOG.md")
    docs_path = Path("docs/markdown/changelog.md")

    # Read and process the changelog
    with open(changelog_path) as f:
        changelog = f.read()

    # Write the combined file
    with open(docs_path, "w") as f:
        f.write(
            f"""\
---
title: MultiQC Version History
description: MultiQC version history and changes
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated. Do not edit the markdown, it will be overwritten.

File path for the source of this content: CHANGELOG.md
~~~~~~~~~~~~~~~~~~~~~~~
-->
        
{changelog}"""
        )

    print(f"Generated {docs_path}")
=======
    print(f"Generated {mdx_path}")
>>>>>>> 314d817d3 (Autogenerate docs into cloned seqeralabs-docs):scripts/make_docs.py
