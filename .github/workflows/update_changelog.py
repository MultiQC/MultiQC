"""
From a PR title and description, add a line into the CHANGELOG.

If a PR starts with "New module: ", add a line under "New modules" in the changelog.

If a PR starts with a name of an existing module, add a line under "Module updates" in the changelog.

Everything else will go under "MultiQC updates" in the changelog, unless it ends with "(chore)" or "(docs)".
"""

import os
import re
import sys
from pathlib import Path

# Get the PR title and description from the environment
pr_title = os.environ["PR_TITLE"]
pr_number = os.environ["PR_NUMBER"]

# Trim the PR number appended to the end of the tile, e.g. "... (#2026)".
if not pr_title.endswith(f" (#{pr_number})"):
    print(f"PR title '{pr_title}' doesn't end with PR number {pr_number}", file=sys.stderr)
    sys.exit(1)
pr_title = pr_title.removesuffix(f" (#{pr_number})")


# Appending "(chore)" or "(docs)" to the PR title will skip logging this change.
if pr_title.endswith("(chore)") or pr_title.endswith("(docs)"):
    section = None
    print("Skipping logging this change as it's a chore or docs update")
    sys.exit(0)


def find_module_info(module_name):
    module_name = module_name.lower()
    modules_dir = Path("../../multiqc/modules")
    py_path = None
    for dir_name in os.listdir(modules_dir):
        if dir_name.lower() == module_name:
            module_dir = modules_dir / dir_name
            py_path = module_dir / f"{dir_name}.py"
            if not py_path.exists():
                print(f"Folder for {module_name} exists, but doesn't have a {py_path} file", file=sys.stderr)
                sys.exit(1)
            break

    if not py_path:  # Module not found
        return None
    with py_path.open("r") as f:
        contents = f.read()
    if not (m := re.search(r'name="([^"]+)"', contents)):
        return None
    name = m.group(1)
    if not (m := re.search(r'href="([^"]+)"', contents)):
        return None
    url = m.group(1)
    if not (m := re.search(r'info="([^"]+)"', contents)):
        if not (m := re.search(r'info="""([^"]+)"""', contents)):
            return None
    info = m.group(1)
    # Strip consecutive spaces and newlines.
    info = re.sub(r"\s+", " ", info)
    return {"name": name, "url": url, "info": info}


mod = None
section = "### MultiQC updates"  # Default section for non-module updates.
if pr_title.lower().startswith("new module: "):
    # PR introduces a new module.
    section = "### New modules"
    module_name = pr_title.split(":")[1].strip()
    mod = find_module_info(module_name)
    if not mod:
        # That should normally never happen because the other CI would fail and block
        # merging of the PR.
        print(
            f"Cannot load a module with name {module_name}",
            file=sys.stderr,
        )
        sys.exit(1)
else:
    # Checking if it's an existing module update.
    maybe_mod_name = pr_title.split(":")[0]
    mod = find_module_info(maybe_mod_name)
    if mod is not None:
        section = "### Module updates"


# Building the change log line for the PR.
pr_link = f"([#{pr_number}](https://github.com/ewels/MultiQC/pull/{pr_number}))"
if section == "### New modules":
    new_lines = [
        f"- [**{mod['name']}**]({mod['url']})\n",
        f"  - {mod['name']} {mod['info']}\n",
    ]
elif section == "### Module updates":
    assert mod is not None
    new_lines = [
        f"- **{mod['name']}**\n",
        f"  - {pr_title} {pr_link}\n",
    ]
else:
    new_lines = [
        f"- {pr_title} {pr_link}\n",
    ]

# Get the current changelog lines. We will print them back as is, except for one new
# line corresponding to this new PR.
with open("../../CHANGELOG.md", "r") as f:
    changelog = f.readlines()
new_changelog = []

# Find the next line in the change log that matches the pattern "## MultiQC v.*dev"
# If it doesn't exist, exist with code 1 (let's assume that a new section is added
# manually or by CI when a release is pushed).
# Else, find the next line that matches the `section` variable, and insert a new line
# under it (we also assume that section headers are added already).
inside_version_dev = False
after_version_dev = False
while changelog:
    line = changelog.pop(0)

    if line.startswith("## "):
        new_changelog.append(line)
        if after_version_dev:
            continue

        # Parse version from the line ## MultiQC v1.10dev or
        # ## [MultiQC v1.15](https://github.com/ewels/MultiQC/releases/tag/v1.15) ...
        m = re.match(r".*MultiQC (v\d+\.\d+(dev)?).*", line)
        if m is None:
            print(f"Cannot parse version from line {line.strip()}.", file=sys.stderr)
            sys.exit(1)
        version = m.group(1)

        if not inside_version_dev and not version.endswith("dev"):
            print(
                "Can't find a 'dev' version section in the changelog. Make sure "
                "it's created, and sections MultiQC updates, New modules and "
                "Module updates are added under it.",
                file=sys.stderr,
            )
            sys.exit(1)
        if inside_version_dev:
            if version.endswith("dev"):
                print(
                    f"Found another 'dev' version section in the changelog, make"
                    f"sure to change it to a 'release' stable version tag. "
                    f"Line: {line.strip()}",
                    file=sys.stderr,
                )
                sys.exit(1)
            inside_version_dev = False
            after_version_dev = True
        else:
            inside_version_dev = True

    elif inside_version_dev and line.lower().startswith(section.lower()):
        if new_lines is None:
            print(f"Already added new lines into section {section}, is the " f"section duplicated?", file=sys.stderr)
            sys.exit(1)
        new_changelog.append(line)
        # Collecting lines until the next section.
        section_lines = []
        while True:
            line = changelog.pop(0)
            if line.startswith("##"):
                new_changelog.append("\n")
                new_changelog.extend(section_lines)
                new_changelog.extend(new_lines)
                new_changelog.append("\n")
                print(f"Updated CHANGELOG.md section '{section}' with lines: {new_lines}")
                new_lines = None
                # pushing back the next section header line
                changelog.insert(0, line)
                break
            elif line.strip():
                section_lines.append(line)
    else:
        new_changelog.append(line)


# Get the current changelog lines. We will print them back as is, except for one new
# line corresponding to this new PR.
with open("../../CHANGELOG.md", "w") as f:
    f.writelines(new_changelog)
