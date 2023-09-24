"""
To be called by a CI action, assumes PR_TITLE and PR_NUMBER, and GITHUB_WORKSPACE environment variables are set.

Adds a line into the CHANGELOG.md:
If a PR title starts with "New module: ", adds a line under the ""### New modules" section.
If a PR starts with a name of an existing module, adds a line under "### Module updates".
Everything else will go under "MultiQC updates" in the changelog, unless "(chore)" or "(docs)" is appended to the title.

Other assumptions:
- CHANGELOG.md has a running section for an ongoing "dev" version (i.e. titled "## MultiQC vX.Ydev").
- Under that section, there are sections "### MultiQC updates", "### New modules" and "### Module updates".
- For module meta info, checks the file multiqc/modules/<module_name>/<module_name>.py.
"""

import os
import re
import sys
from pathlib import Path

REPO_URL = "https://github.com/ewels/MultiQC"

# Assumes the environment is set by the GitHub action.
pr_title = os.environ["PR_TITLE"]
pr_number = os.environ["PR_NUMBER"]
comment = os.environ.get("COMMENT", "")
base_path = Path(os.environ.get("GITHUB_WORKSPACE", ""))

# Trim the PR number added when GitHub squashes commits, e.g. "Module: Updated (#2026)"
pr_title = pr_title.removesuffix(f" (#{pr_number})")

changelog_path = base_path / "CHANGELOG.md"


def find_module_info(module_name):
    """
    Helper function to load module meta info. With current setup, can't really just
    import the module and call `mod.info`, as the module does the heavy work on
    initialization. But that's actually alright: we avoid installing and importing
    MultiQC and the action runs faster.
    """
    module_name = module_name.lower()
    modules_dir = base_path / "multiqc/modules"
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
    # Reduce consecutive spaces and newlines.
    info = re.sub(r"\s+", " ", info)
    return {"name": name, "url": url, "info": info}


# Determine the type of the PR: new module, module update, or core update.
mod = None
section = "### MultiQC updates"  # Default section for non-module (core) updates.
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
        pr_title = pr_title.split(":")[1].strip().capitalize()

# Now that we determined the PR type, preparing the change log entry.
pr_link = f"([#{pr_number}]({REPO_URL}/pull/{pr_number}))"
if comment := comment.removeprefix("/changelog ").strip():
    new_lines = [
        f"- {comment} {pr_link}\n",
    ]
elif section == "### New modules":
    new_lines = [
        f"- [**{mod['name']}**]({mod['url']}) {pr_link}\n",
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

# Finally, updating the changelog.
# Read the current changelog lines. We will print them back as is, except for one new
# entry, corresponding to this new PR.
with changelog_path.open("r") as f:
    orig_lines = f.readlines()
updated_lines = []

# Find the next line in the change log that matches the pattern "## MultiQC v.*dev"
# If it doesn't exist, exist with code 1 (let's assume that a new section is added
# manually or by CI when a release is pushed).
# Else, find the next line that matches the `section` variable, and insert a new line
# under it (we also assume that section headers are added already).
inside_version_dev = False
while orig_lines:
    line = orig_lines.pop(0)

    if line.startswith("## "):  # Version header, e.g. "## MultiQC v1.10dev"
        updated_lines.append(line)

        # Parse version from the line ## MultiQC v1.10dev or
        # ## [MultiQC v1.15](https://github.com/ewels/MultiQC/releases/tag/v1.15) ...
        if not (m := re.match(r".*MultiQC (v\d+\.\d+(dev)?).*", line)):
            print(f"Cannot parse version from line {line.strip()}.", file=sys.stderr)
            sys.exit(1)
        version = m.group(1)

        if not inside_version_dev:
            if not version.endswith("dev"):
                print(
                    "Can't find a 'dev' version section in the changelog. Make sure "
                    "it's created, and sections MultiQC updates, New modules and "
                    "Module updates are added under it.",
                    file=sys.stderr,
                )
                sys.exit(1)
            inside_version_dev = True
        else:
            if version.endswith("dev"):
                print(
                    f"Found another 'dev' version section in the changelog, make"
                    f"sure to change it to a 'release' stable version tag. "
                    f"Line: {line.strip()}",
                    file=sys.stderr,
                )
                sys.exit(1)
            # We are past the dev version, so just add back the rest of the lines and break.
            updated_lines.extend(orig_lines)
            break

    elif inside_version_dev and line.lower().startswith(section.lower()):  # Section of interest header
        if new_lines is None:
            print(f"Already added new lines into section {section}, is the section duplicated?", file=sys.stderr)
            sys.exit(1)
        updated_lines.append(line)
        # Collecting lines until the next section.
        section_lines = []
        while True:
            line = orig_lines.pop(0)
            if line.startswith("##"):
                # Found the next section header, so need to put all the lines we collected.
                updated_lines.append("\n")
                updated_lines.extend(section_lines)
                updated_lines.extend(new_lines)
                updated_lines.append("\n")
                print(f"Updated {changelog_path} section '{section}' with lines:\n" + "".join(new_lines))
                new_lines = None
                # Pushing back the next section header line
                orig_lines.insert(0, line)
                break
            elif line.strip():
                # if the line already contains a link to the PR, don't add it again.
                if line.strip().endswith(pr_link):
                    existing = line + "".join(orig_lines[: len(new_lines) - 1])
                    if "".join(new_lines) == existing:
                        print(f"Found existing identical entry for this pull request #{pr_number}:")
                        print(existing)
                        sys.exit(0)
                    else:
                        print(f"Found existing entry for this pull request #{pr_number}. It will be replaced:")
                        print(existing)
                        for _ in range(len(new_lines) - 1):
                            orig_lines.pop(0)
                else:
                    section_lines.append(line)
    else:
        updated_lines.append(line)


# Finally, writing the updated lines back.
with changelog_path.open("w") as f:
    f.writelines(updated_lines)
