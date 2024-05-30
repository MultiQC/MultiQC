"""
Script to prepare a release changelog section.

Takes the git commit history since the last tag, assuming each commit message is a changelog
entry (if it's not suffixed with [no changelog]), and creates a markdown.
"""

import datetime
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import importlib_metadata
import yaml

REPO_URL = "https://github.com/MultiQC/MultiQC"
WORKSPACE_PATH = Path(os.environ.get("GITHUB_WORKSPACE", "."))
MODULES_SUBDIR = Path("multiqc/modules")


def main():
    entries_by_section = {
        "### MultiQC updates": [],
        "### New modules": [],
        "### Module updates": defaultdict(list),
    }

    last_tag = run_cmd(f"cd {WORKSPACE_PATH} && git describe --tags --abbrev=0").stdout.strip()
    commit_history = run_cmd(f"cd {WORKSPACE_PATH} && git log --oneline --decorate='' {last_tag}..HEAD").stdout
    for commit in commit_history.splitlines():
        sha, message = commit.split(" ", 1)
        if skip_commit(message):
            continue
        if not (m := re.search(r"(.+) \(#(\d+)\)$", message)):
            continue
        message = m.group(1)
        pr_number = m.group(2)
        pr_link = f"([#{pr_number}]({REPO_URL}/pull/{pr_number}))"

        section, mod_info = determine_change_type(message, sha)
        if section == "### New modules":
            entries_by_section[section].append((mod_info, message, pr_link))
        elif section == "### Module updates":
            assert mod_info is not None
            if ":" in message:
                message = message.split(":", maxsplit=1)[1].strip()
            entries_by_section[section][mod_info["name"]].append((message, pr_link))
        else:  # "### MultiQC updates"
            entries_by_section[section].append((message, pr_link))
    run_cmd(f"cd {WORKSPACE_PATH} && git checkout main")

    dev_version = importlib_metadata.version("multiqc")  # 1.22.dev0
    new_version = dev_version.removesuffix(".dev0")
    print(f"## [MultiQC v{new_version}]({REPO_URL}/releases/tag/v{new_version}) - {datetime.date.today().isoformat()}")
    print("")
    if entries_by_section["### MultiQC updates"]:
        print("### MultiQC updates")
        print("")
        for entry in entries_by_section["### MultiQC updates"]:
            message, pr_link = entry
            print(f"- {message} {pr_link}")
        print("")
    if entries_by_section["### New modules"]:
        print("### New modules")
        print("")
        for entry in entries_by_section["### New modules"]:
            mod_info, descr, pr_link = entry
            print(f"- [**{mod_info['name']}**]({mod_info['url']}) {pr_link}")
            print(f"  - {descr} {mod_info['info']}")
        print("")
    if entries_by_section["### Module updates"]:
        print("### Module updates")
        print("")
        for mod_name, entries in entries_by_section["### Module updates"].items():
            print("- **" + mod_name + "**")
            for entry in entries:
                descr, pr_link = entry
                print(f"  - {descr} {pr_link}")
        print("")


def skip_commit(message: str) -> bool:
    return any(
        line in message.lower()
        for line in [
            "skip changelog",
            "skip change log",
            "no changelog",
            "no change log",
            "bump version",
        ]
    )


def run_cmd(cmd):
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error executing command: {result.stderr}")
    return result


def diff_for_a_file(commit_sha, path) -> str:
    """
    Returns the diff for a specific file altered in the PR.
    """
    result = run_cmd(f"cd {WORKSPACE_PATH} && git show {commit_sha}")

    lines = result.stdout.splitlines()
    while lines:
        line = lines.pop(0)
        if line.startswith("index "):
            a = lines.pop(0)
            b = lines.pop(0)
            if a == f"--- a/{path}" and b == f"+++ b/{path}":
                diff = []
                while lines:
                    line = lines.pop(0)
                    if line.startswith("diff "):
                        break
                    diff.append(line)
                return "\n".join(diff)


def determine_change_type(message, commit_sha) -> tuple[str, dict]:
    """
    Determine the type of the PR: new module, module update, or core update.
    Returns a tuple of the section name and the module info.
    """
    added_modules = modules_added_by_commit(commit_sha)
    modified_modules = modules_modified_by_commit(commit_sha)
    modified_modules = {
        m for m in modified_modules if m not in {"profile_runtime.py", "custom_content.py", "base_module.py"}
    }

    # Sanity check PR name suggesting a newly added module
    if message.lower().capitalize().startswith("New module: "):
        if len(added_modules) == 0:
            raise RuntimeError(
                f"Could not find a new folder in '{MODULES_SUBDIR}' with expected python files for the new module"
            )
        if len(added_modules) > 1:
            RuntimeError(f"Found multiple added modules: {added_modules}")

    if len(added_modules) == 1 and len(modified_modules) == 0:
        mod_info = find_module_info(added_modules[0])
        proper_pr_title = f"New module: {mod_info['name']}"
        if message != proper_pr_title:
            try:
                run_cmd(f"cd {WORKSPACE_PATH} && gh pr edit --title '{proper_pr_title}'")
            except (RuntimeError, subprocess.CalledProcessError) as e:
                print(
                    f"Error executing command: {e}. Please alter the title manually: '{proper_pr_title}'",
                    file=sys.stderr,
                )
        return "### New modules", mod_info

    elif len(modified_modules) == 1 and len(added_modules) == 0:
        mod_info = find_module_info(modified_modules.pop())
        return "### Module updates", mod_info

    section = "### MultiQC updates"  # Default section for non-module (core) updates.
    return section, {}


def find_module_info(mod_name: str) -> dict[str]:
    """
    Helper function to load module meta info. With current setup, can't really just
    import the module and call `mod.info`, as the module does the heavy work on
    initialization. But that's actually alright: we avoid installing and importing
    MultiQC and the action runs faster.
    """
    py_path = WORKSPACE_PATH / MODULES_SUBDIR / mod_name / f"{mod_name}.py"
    if not py_path.exists():
        return {}

    if py_path.name == "custom_content.py":
        return {"name": "Custom content", "anchor": "custom-content", "url": "", "info": ""}

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


def files_altered_by_pr(commit_sha, types=None) -> set[Path]:
    """
    Returns a list of files added or modified by the PR (depending on `types`,
    which can be a subset of `{'added', 'modified'}`)
    """
    if types is None:
        types = {"added"}

    result = run_cmd(f"cd {WORKSPACE_PATH} && git show {commit_sha}")

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


def modules_added_by_commit(commit_sha) -> list[str]:
    """
    Returns paths to the modules added by the PR.
    """
    added_modules = []
    altered_files = files_altered_by_pr(commit_sha, {"added"})
    for path in altered_files:
        if path.name == "__init__.py" and str(path).startswith(f"{MODULES_SUBDIR}/"):
            mod_name = path.parent.name
            if (mod_py := path.parent / f"{mod_name}.py") in altered_files:
                if (WORKSPACE_PATH / mod_py).exists():
                    added_modules.append(mod_name)
    return added_modules


def modules_modified_by_commit(commit_sha) -> set[str]:
    """
    Returns paths to the "<module>.py" files of the altered modules.
    """
    altered_files = files_altered_by_pr(commit_sha, {"modified"})

    # First, special case for search patterns.
    keys_added_in_search_patterns = set()
    keys_modified_in_search_patterns = set()
    sp_paths = [f for f in altered_files if f.name == "search_patterns.yaml"]
    if sp_paths:
        sp_path = sp_paths[0]
        # Find modules that changed, e.g. collect "htseq":
        # htseq:
        #   - contents_re: '^(feature\tcount|\w+\t\d+)$'
        #   + contents_re: '^(feature\tcount|\w+.*\t\d+)$'
        #   num_lines: 1
        # Checkout code before commit
        run_cmd(f"cd {WORKSPACE_PATH} && git checkout {commit_sha}^1")
        with (WORKSPACE_PATH / sp_path).open() as f:
            old_text = f.read()
        # Get contents of the file changed by the PR.
        run_cmd(f"cd {WORKSPACE_PATH} && git checkout {commit_sha}")
        with (WORKSPACE_PATH / sp_path).open() as f:
            new_text = f.read()
        if old_text != new_text:
            old_data = yaml.safe_load(old_text)
            new_data = yaml.safe_load(new_text)
            for new_skey in new_data:
                # Added module?
                if new_skey not in old_data:
                    keys_added_in_search_patterns.add(new_skey)
            for old_skey in old_data:
                # Removed module?
                if old_skey not in new_data:
                    keys_modified_in_search_patterns.add(old_skey)
                # Modified module?
                elif old_data.get(old_skey) != new_data.get(old_skey):
                    keys_modified_in_search_patterns.add(old_skey)

    mod_names = set()
    if keys_modified_in_search_patterns:
        mod_name = keys_modified_in_search_patterns.pop().split("/")[0]
        mod_names.add(mod_name)

    # Now adding module-specific files
    for path in altered_files:
        if path.is_relative_to(MODULES_SUBDIR):
            mod_name = path.relative_to(MODULES_SUBDIR).parts[0]
            mod_names.add(mod_name)

    return mod_names


if __name__ == "__main__":
    main()
