"""
To be called by a CI action. Assumes the following environment variables are set:
PR_TITLE, PR_NUMBER, GITHUB_WORKSPACE.

Adds a line into the CHANGELOG.md:
* If a PR title starts with "New module: ", checks that a single module is added,
 and appends an entry under the ""### New modules" section. 
* If a single module was modified, checks that the PR starts with a name of the
 modified module (e.g. "FastQC: new stuff") and adds a line under "### Module updates".
* All other change will go under the "### MultiQC updates" section.
* If an entry for the PR is already added, will replace it.

Other assumptions:
- CHANGELOG.md has a running section for an ongoing "dev" version 
(i.e. titled "## MultiQC vX.Ydev").
- Under that section, there are sections "### MultiQC updates", "### New modules" 
and "### Module updates".
- For module's info, checks the file multiqc/modules/<module_name>/<module_name>.py.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

REPO_URL = "https://github.com/ewels/MultiQC"
MODULES_DIR = "multiqc/modules"

# Assumes the environment is set by the GitHub action.
pr_title = os.environ["PR_TITLE"]
pr_number = os.environ["PR_NUMBER"]
comment = os.environ.get("COMMENT", "")
workspace_path = Path(os.environ.get("GITHUB_WORKSPACE", ""))

assert pr_title, pr_title
assert pr_number, pr_number

# Trim the PR number added when GitHub squashes commits, e.g. "Module: Updated (#2026)"
pr_title = pr_title.removesuffix(f" (#{pr_number})")

changelog_path = workspace_path / "CHANGELOG.md"


def _find_module_info(py_path: Path) -> dict[str]:
    """
    Helper function to load module meta info. With current setup, can't really just
    import the module and call `mod.info`, as the module does the heavy work on
    initialization. But that's actually alright: we avoid installing and importing
    MultiQC and the action runs faster.
    """
    with py_path.open("r") as f:
        contents = f.read()

    if not (m := re.search(r'name="([^"]+)"', contents)):
        return {}
    name = m.group(1)

    if not (m := re.search(r'anchor="([^"]+)"', contents)):
        return {}
    anchor = m.group(1)

    if not (m := re.search(r'href="([^"]+)"', contents)):
        return {}
    url = m.group(1)

    if not (m := re.search(r'info="([^"]+)"', contents)):
        if not (m := re.search(r'info="""([^"]+)"""', contents)):
            return {}
    info = m.group(1)

    # Reduce consecutive spaces and newlines.
    info = re.sub(r"\s+", " ", info)
    return {"name": name, "anchor": anchor, "url": url, "info": info}


def _files_altered_by_pr(pr_number, types=None) -> set[Path]:
    """
    Returns a list of files added or modified by the PR (depending on `types`,
    which can be a subset of `{'added', 'modified'}`)
    """
    if types is None:
        types = {"added"}

    cmd = f"cd {workspace_path} && gh pr diff {pr_number}"
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error executing command: {result.stderr}")

    paths = set()
    lines = result.stdout.splitlines()
    while lines:
        line = lines.pop(0)
        if line.startswith("index "):
            a = lines.pop(0)
            b = lines.pop(0)
            if a.startswith("--- /dev/null") and b.startswith("+++ b/") and "added" in types:
                paths.add(Path(b.removeprefix("+++ b/")))
            elif a.startswith("--- a/") and b.startswith("+++ /dev/null") and "deleted" in types:
                paths.add(Path(a.removeprefix("--- a/")))
            elif a.startswith("--- a/") and b.startswith("+++ b/") and "modified" in types:
                paths.add(Path(a.removeprefix("--- a/")))
    return paths


def _modules_added_by_pr(pr_number) -> list[Path]:
    """
    Returns paths to the modules added by the PR.
    """
    mod_py_files = []
    altered_files = _files_altered_by_pr(pr_number, {"added"})
    for path in altered_files:
        if path.name == "__init__.py" and str(path).startswith(f"{MODULES_DIR}/"):
            mod_anchor = path.parent.name
            if (mod_py := path.parent / f"{mod_anchor}.py") in altered_files:
                if (mod_py := workspace_path / mod_py).exists():
                    mod_py_files.append(mod_py)
    return mod_py_files


def _modules_modified_by_pr(pr_number) -> list[Path]:
    """
    Returns paths to the modules modified by the PR.
    """
    mod_py_files = []
    altered_files = _files_altered_by_pr(pr_number, {"modified"})
    for path in altered_files:
        if str(path).startswith(f"{MODULES_DIR}/"):
            mod_anchor = path.parent.name
            mod_py = path.parent / f"{mod_anchor}.py"
            if (mod_py := workspace_path / mod_py).exists():
                mod_py_files.append(mod_py)
    return mod_py_files


def _determine_change_type(pr_title, pr_number) -> tuple[str, dict]:
    """
    Determine the type of the PR: new module, module update, or core update.
    Returns a tuple of the section name and the module info.
    """

    if pr_title.lower().capitalize().startswith("New module: "):
        mod_py_files = _modules_added_by_pr(pr_number)
        if len(mod_py_files) == 0:
            raise RuntimeError(
                f"Could not find a new folder in '{MODULES_DIR}' with expected python files for the new module"
            )
        if len(mod_py_files) > 1:
            RuntimeError(f"Found multiple added modules: {mod_py_files}")
        else:
            mod_info = _find_module_info(mod_py_files[0])
            proper_pr_title = f"New module: {mod_info['name']}"
            if pr_title != proper_pr_title:
                cmd = f"cd {workspace_path}; gh pr edit --title '{proper_pr_title}'"
                print(cmd)
                try:
                    subprocess.run(cmd, shell=True)
                except subprocess.CalledProcessError as e:
                    print(
                        f"Error executing command: {e}. Please alter the title manually: '{proper_pr_title}'",
                        file=sys.stderr,
                    )
            return "### New modules", mod_info

    # Check what modules were changed by the PR, and if the title starts with on
    # of the names of the changed modules, assume it's a module update.
    modified_mod_py_files = _modules_modified_by_pr(pr_number)
    if len(modified_mod_py_files) == 1:
        mod_info = _find_module_info(modified_mod_py_files.pop())
        if pr_title.lower().startswith(f"{mod_info['name'].lower()}: "):
            return "### Module updates", mod_info

    section = "### MultiQC updates"  # Default section for non-module (core) updates.
    return section, {}


# Determine the type of the PR: new module, module update, or core update.
section, mod = _determine_change_type(pr_title, pr_number)

# Prepare the change log entry.
pr_link = f"([#{pr_number}]({REPO_URL}/pull/{pr_number}))"
if comment := comment.removeprefix("@multiqc-bot changelog").strip().capitalize():
    new_lines = [
        f"- {comment.capitalize()} {pr_link}\n",
    ]
elif section == "### New modules":
    new_lines = [
        f"- [**{mod['name']}**]({mod['url']}) {pr_link}\n",
        f"  - {mod['name']} {mod['info']}\n",
    ]
elif section == "### Module updates":
    assert mod is not None
    descr = pr_title.split(":", maxsplit=1)[1].strip().capitalize()
    new_lines = [
        f"- **{mod['name']}**: {descr} {pr_link}\n",
    ]
else:
    new_lines = [
        f"- {pr_title.capitalize()} {pr_link}\n",
    ]


# Finally, updating the changelog.
# Read the current changelog lines. We will print them back as is, except for one new
# entry, corresponding to this new PR.
with changelog_path.open("r") as f:
    orig_lines = f.readlines()
updated_lines = []


def _skip_existing_entry_for_this_pr(line, same_section=True):
    if line.strip().endswith(pr_link):
        existing = line + "".join(orig_lines[: len(new_lines) - 1])
        if "".join(new_lines) == existing and same_section:
            print(f"Found existing identical entry for this pull request #{pr_number}:")
            print(existing)
            sys.exit(0)
        else:
            print(
                f"Found existing entry for this pull request #{pr_number}. It will be replaced and/or moved to proper section"
            )
            print(existing)
            for _ in range(len(new_lines)):
                try:
                    line = orig_lines.pop(0)
                except IndexError:
                    break
    return line


# Find the next line in the change log that matches the pattern "## MultiQC v.*dev"
# If it doesn't exist, exist with code 1 (let's assume that a new section is added
# manually or by CI when a release is pushed).
# Else, find the next line that matches the `section` variable, and insert a new line
# under it (we also assume that section headers are added already).
inside_version_dev = False
already_added_entry = False
while orig_lines:
    line = orig_lines.pop(0)

    # If the line already contains a link to the PR, don't add it again.
    line = _skip_existing_entry_for_this_pr(line, same_section=False)

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
            while orig_lines:
                line = orig_lines.pop(0)
                line = _skip_existing_entry_for_this_pr(line, same_section=False)
                if line:
                    updated_lines.append(line)
            break
        continue

    if inside_version_dev and line.lower().startswith(section.lower()):  # Section of interest header
        if already_added_entry:
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
                _updated_lines = [_l for _l in section_lines + new_lines if _l.strip()]
                if section == "### Module updates":
                    _updated_lines = sorted(_updated_lines)
                updated_lines.extend(_updated_lines)
                updated_lines.append("\n")
                print(f"Updated {changelog_path} section '{section}' with lines:\n" + "".join(new_lines))
                already_added_entry = True
                # Pushing back the next section header line
                orig_lines.insert(0, line)
                break
            # If the line already contains a link to the PR, don't add it again.
            line = _skip_existing_entry_for_this_pr(line, same_section=True)
            section_lines.append(line)
    else:
        updated_lines.append(line)


# Finally, writing the updated lines back.
with changelog_path.open("w") as f:
    f.writelines(updated_lines)
