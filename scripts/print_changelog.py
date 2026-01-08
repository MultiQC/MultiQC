#!/usr/bin/env python

"""
Script to prepare a release changelog section.

Takes the GitHub PR history since the last milestone and prints a markdown.
"""

import datetime
import os
import subprocess
from pathlib import Path
from typing import Dict, List

from github import Github
from github.PullRequest import PullRequest

REPO_ID = "MultiQC/MultiQC"
REPO_URL = f"https://github.com/{REPO_ID}"
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")
WORKSPACE_PATH = Path(os.environ.get("GITHUB_WORKSPACE", "."))
MODULES_SUBDIR = Path("multiqc/modules")

if not GITHUB_TOKEN:
    raise ValueError("Please set the GITHUB_TOKEN environment variable")


def run_cmd(cmd):
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error executing command: {result.stderr}")
    return result


def get_milestone_prs(repo, current_tag: str, limit=100) -> List[PullRequest]:
    """
    Get PRs for the current milestone using direct milestone filtering.

    This version directly queries PRs by milestone, ensuring completeness
    and accuracy while minimizing API calls by using the Issue API.
    """
    # Find the milestone object for the current tag
    milestones = {m.title: m for m in repo.get_milestones(state="all")}

    current_milestone = milestones.get(current_tag)
    if not current_milestone:
        raise ValueError(f"Current milestone '{current_tag}' not found")

    print(
        f"Found milestone '{current_tag}' with {current_milestone.open_issues + current_milestone.closed_issues} total issues/PRs"
    )

    # Get issues/PRs for this specific milestone (GitHub API treats PRs as issues)
    # This is more reliable than iterating through all PRs
    issues = repo.get_issues(state="closed", milestone=current_milestone, sort="updated", direction="desc")

    # Filter to merged PRs only
    all_pulls = []
    for issue in issues:
        if issue.pull_request is not None:  # This is a PR
            try:
                pr = repo.get_pull(issue.number)
                if pr.merged:  # Only include merged PRs
                    all_pulls.append(pr)
                    if len(all_pulls) >= limit:
                        print(f"Reached limit of {limit} PRs")
                        break
            except Exception as e:
                print(f"Error fetching PR #{issue.number}: {e}")
                continue

    print(f"Found {len(all_pulls)} merged PRs for milestone '{current_tag}'")
    return all_pulls


def assert_milestone_exists(milestones, tag: str):
    for m in milestones:
        if m.title == tag:
            return True
    raise ValueError(
        f"Could not find milestone '{tag}'. Please double check that you updated the "
        f"version tag in pyproject.toml and the milestone name in the GitHub repo matches it."
    )


def get_version_from_tag() -> str:
    with open("pyproject.toml") as f:
        for line in f:
            if line.startswith("version ="):
                version = line.split(" = ")[1].strip().strip('"')
                return version.removesuffix("dev")
    raise ValueError("Could not find version in pyproject.toml")


def main():
    current_tag = f"v{get_version_from_tag()}"
    previous_minor_tag = run_cmd(f"cd {WORKSPACE_PATH} && git describe --tags --abbrev=0").stdout.strip()
    if previous_minor_tag.count(".") == 2:  # 1.24.1 -> 1.24
        previous_minor_tag = previous_minor_tag.rsplit(".", 1)[0]
    repo = Github(login_or_token=GITHUB_TOKEN).get_repo(REPO_ID)
    milestones = repo.get_milestones(state="all")
    assert_milestone_exists(milestones, current_tag)
    assert_milestone_exists(milestones, previous_minor_tag)
    prs: List[PullRequest] = get_milestone_prs(repo, current_tag)

    label_to_section: Dict[str, str] = {
        "module: new": "New modules",
        "bug: module": "Module fixes",
        "module: enhancement": "Module updates",
        "module: change": "Module updates",
        "bug: core": "Fixes",
        "core: infrastructure": "Infrastructure and packaging",
        "core: refactoring": "Optimization, refactoring and typing",
        "documentation": "Chores",
    }
    sections_to_prs: Dict[str, List[PullRequest]] = {
        "New modules": [],
        "Feature updates and improvements": [],
        "Module updates": [],
        "Fixes": [],
        "Module fixes": [],
        "Refactoring and typing": [],
        "Infrastructure and packaging": [],
        "Optimization, refactoring and typing": [],
        "Chores": [],
    }
    for pr in prs:
        if skip_pr(pr.title):
            continue
        if pr.labels:
            for label in pr.labels:
                section_name = label_to_section.get(label.name, "Feature updates and improvements")
                sections_to_prs[section_name].append(pr)
        else:
            sections_to_prs["Feature updates and improvements"].append(pr)

    print("")
    print(f"## [MultiQC {current_tag}]({REPO_URL}/releases/tag/{current_tag}) - {datetime.date.today().isoformat()}")
    print("")

    for section, prs in sections_to_prs.items():
        if section == "Chores":
            continue
        print(f"### {section}")
        print("")
        for pr in prs:
            link = f"([#{pr.number}]({REPO_URL}/pull/{pr.number}))"
            print(f"- {pr.title} {link}")
        print("")


def skip_pr(message: str) -> bool:
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


if __name__ == "__main__":
    main()
