import re
from pathlib import Path
import subprocess
from typing import TypedDict


def extract_first_second_level_section(markdown) -> str:
    # Match the first second-level section and capture its content
    pattern = r"(##\s+.*?)(?=(\n##\s+|$))"
    match = re.search(pattern, markdown, re.DOTALL)

    if not match:
        raise ValueError("Could not find first second-level section in changelog")

    # Return the matched section
    return match.group(1).strip()


class ChangelogData(TypedDict):
    version: str
    url: str
    date: str
    summary: str
    the_rest: str


def extract_latest_version(changelog_content) -> ChangelogData:
    """Extract the latest version section from the changelog."""
    # Skip the first header line
    last_version_changelog = extract_first_second_level_section(changelog_content)

    # Find the first version section and extract the version, url and date
    matches = re.search(r"^## \[MultiQC v([\d.]+)\]\((.*?)\) - ([\d-]+)", last_version_changelog)
    if not matches:
        raise ValueError("Could not find version information in changelog")
    version, url, date = matches.groups()

    # Remove the first line
    last_version_changelog = "\n".join(last_version_changelog.splitlines()[1:]).strip()
    # Get data until the third-level header
    next_header_index = last_version_changelog.find("###")
    if next_header_index == -1:
        raise ValueError("Could not find next header in changelog")
    summary = last_version_changelog[:next_header_index].strip()
    the_rest = last_version_changelog[next_header_index:].strip()

    return {"version": version, "date": date, "url": url, "summary": summary, "the_rest": the_rest}


# Read CHANGELOG.md
changelog_path = Path("CHANGELOG.md")
if not changelog_path.exists():
    raise FileNotFoundError("CHANGELOG.md not found")

# Create seqeralabs-docs directory if it doesn't exist
docs_dir = Path("seqeralabs-docs")
# Clone seqeralabs-docs repo if it doesn't exist
if not docs_dir.exists():
    print("Cloning seqeralabs-docs repository...")
    repo_url = "https://github.com/seqeralabs/docs.git"
    try:
        subprocess.run(["git", "clone", repo_url, str(docs_dir)], check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to clone repository: {e}")

# Extract latest version
with open(changelog_path) as f:
    changelog_data: ChangelogData = extract_latest_version(f.read())

# Create output directory
output_dir = docs_dir / "changelog" / "multiqc"
output_dir.mkdir(parents=True, exist_ok=True)

# Create output file
mdx_content: str = f"""---
title: MultiQC v{changelog_data["version"]}
date: {changelog_data["date"]}
tags: [multiqc]
---

{changelog_data["summary"]}

{{/* truncate */}}

{changelog_data["the_rest"]}
"""
with open(output_path := output_dir / f"v{changelog_data['version']}.mdx", "w") as f:
    f.write(mdx_content)

# Print version for GitHub Actions
print(f"::set-output name=version::{changelog_data['version']}")
