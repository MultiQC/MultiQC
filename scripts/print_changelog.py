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
GITHUB_TOKEN = os.environ["GITHUB_TOKEN"]
WORKSPACE_PATH = Path(os.environ.get("GITHUB_WORKSPACE", "."))
MODULES_SUBDIR = Path("multiqc/modules")


def run_cmd(cmd):
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error executing command: {result.stderr}")
    return result


def get_milestone_prs(repo, current_tag: str, previous_tag: str, limit=100) -> List[PullRequest]:
    all_pulls: List[PullRequest] = []

    page = repo.get_pulls(state="closed", sort="updated", direction="desc")
    for p in page:
        if not p.milestone:
            print(f"PR does not have a milestone: {p.number} {p.title}")
        elif p.milestone.title == previous_tag:
            return all_pulls
        elif p.milestone.title == current_tag:
            all_pulls.append(p)
            if len(all_pulls) >= limit:
                print(f"Reached limit of {limit} PRs")
                return all_pulls
        else:
            print(f"PR not in milestone '{p.milestone.title}': {p.number} {p.title}")

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
                return version
    raise ValueError("Could not find version in pyproject.toml")


def main():
    current_tag = f"v{get_version_from_tag()}"
    previous_tag = run_cmd(f"cd {WORKSPACE_PATH} && git describe --tags --abbrev=0").stdout.strip()
    repo = Github(login_or_token=GITHUB_TOKEN).get_repo(REPO_ID)
    milestones = repo.get_milestones(state="all")
    assert_milestone_exists(milestones, current_tag)
    assert_milestone_exists(milestones, previous_tag)
    prs: List[PullRequest] = get_milestone_prs(repo, current_tag, previous_tag)

    label_to_section: Dict[str, str] = {
        "module: new": "New modules",
        "bug: module": "Module fixes",
        "module: enhancement": "Module updates",
        "module: change": "Updates",
        "bug: core": "Fixes",
        "core: infrastructure": "Infrastructure",
        "core: refactoring": "Refactoring",
        "documentation": "Chores",
    }
    sections_to_prs: Dict[str, List[PullRequest]] = {
        "Fixes": [],
        "Updates": [],
        "New modules": [],
        "Module fixes": [],
        "Module updates": [],
        "Refactoring": [],
        "Infrastructure": [],
        "Chores": [],
    }
    for pr in prs:
        if skip_pr(pr.title):
            continue
        if pr.labels:
            for label in pr.labels:
                section_name = label_to_section.get(label.name, "Updates")
                sections_to_prs[section_name].append(pr)
        else:
            sections_to_prs["Updates"].append(pr)

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
