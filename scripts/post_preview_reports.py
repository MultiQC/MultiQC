#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["PyGithub==2.10.0", "requests"]
# ///
"""
Post MultiQC preview report comments on PRs, from builds made by module-report-build.yml.

Run by .github/workflows/module-report.yml. Never runs PR code: it only reads build
results and artifacts through the GitHub API, and comments with the bot token.
"""

import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import requests
from github import Auth, Github
from github.PullRequest import PullRequest
from github.WorkflowRun import WorkflowRun

REPO = os.environ["REPO"]
EVENT = os.environ["EVENT"]
IS_REPORT = bool(re.search(r"(?m)^/report", os.environ.get("COMMENT_BODY", "")))

# Every bot comment contains this, so it can be found and hidden later
MARKER = "MultiQC preview report"
SHA_RE = re.compile(r"<!-- preview-sha: ([0-9a-f]+) -->")
# Keep in sync with module-report-build.yml
MODULE_RE = re.compile(r"^multiqc/modules/([A-Za-z0-9_]+)/")

HOW = """

<details><summary>How this works</summary>

The bot posts a preview report on PRs that edit modules:

1. When the first build finishes, a few minutes after the PR opens. If you push more fixes before then, the report uses the latest commit.
2. Each night, if there are new commits since the last report.
3. When anyone comments `/report`. If the build is still running, the report posts when it finishes.

The bot only posts a new report when there are new commits. If there are none, it reacts 😕 to your comment. The bot hides older report comments.

</details>"""

# Actions token: reads runs and artifacts. Bot PAT: comments, reactions and --attach.
repo = Github(auth=Auth.Token(os.environ["GH_TOKEN"])).get_repo(REPO)
bot = Github(auth=Auth.Token(os.environ["BOT_TOKEN"]))
bot_repo = bot.get_repo(REPO)
BOT_LOGIN = bot.get_user().login


def react(content: str) -> None:
    if IS_REPORT:
        bot_repo.get_issue(int(os.environ["ISSUE"])).get_comment(int(os.environ["COMMENT_ID"])).create_reaction(content)


def post(pr: PullRequest, body: str, attach: Optional[Path] = None) -> None:
    # gh, because the REST API can't upload images
    args = ["gh", "pr", "comment", str(pr.number), "--repo", REPO, "--body", body + HOW]
    if attach:
        args += ["--attach", f"./{attach}"]
    subprocess.run(args, env=dict(os.environ, GH_TOKEN=os.environ["BOT_TOKEN"]), check=True)


def get_targets() -> Tuple[List[PullRequest], str, Optional[WorkflowRun]]:
    """
    Return the PRs to consider, the posting rule, and the build if known.

    Rule "first": post only if the PR has no report yet, or has an unanswered /report.
    Rule "new": post only if the head moved since the last report.
    """
    if EVENT == "issue_comment":
        return [bot_repo.get_pull(int(os.environ["ISSUE"]))], "new", None
    if EVENT == "schedule":
        return list(bot_repo.get_pulls(state="open")), "new", None

    # workflow_run. Fork runs have an empty pull_requests payload, so the build's
    # run-name carries the PR number.
    build = repo.get_workflow_run(int(os.environ["WR_ID"]))
    match = re.match(r"#(\d+) ", build.display_title)
    if not match:
        sys.exit(f"Can't parse PR from run name: {build.display_title}")
    pr = bot_repo.get_pull(int(match.group(1)))
    if build.event != "pull_request":
        return [pr], "new", build
    # Also rejects a fork run naming someone else's PR
    if pr.head.sha != build.head_sha:
        print(f"PR #{pr.number} head moved on or doesn't match this build")
        return [], "first", build
    return [pr], "first", build


def find_build(sha: str) -> Optional[WorkflowRun]:
    runs = repo.get_workflow("module-report-build.yml").get_runs(head_sha=sha, status="completed")
    return next((r for r in runs if r.conclusion in ("success", "failure")), None)


def process(pr: PullRequest, rule: str, build: Optional[WorkflowRun]) -> None:
    sha = pr.head.sha
    modules = sorted({m.group(1) for f in pr.get_files() if (m := MODULE_RE.match(f.filename))})
    if not modules:
        if IS_REPORT:
            post(pr, f"ℹ️ No edited modules found, so no {MARKER}.")
        return

    comments = list(pr.get_issue_comments())
    mine = [c for c in comments if c.user.login == BOT_LOGIN and MARKER in c.body]
    shas = [m.group(1) for c in mine if (m := SHA_RE.search(c.body))]
    last_sha = shas[-1] if shas else ""
    # A /report that arrived while the build was still running, not yet answered by a report
    after_last = comments[comments.index(mine[-1]) + 1 :] if mine else comments
    pending = any(re.search(r"(?m)^/report", c.body) for c in after_last)

    if last_sha == sha:
        react("confused")
        return
    if rule == "first" and last_sha and not pending:
        return

    # No finished build yet: its workflow_run posts the report when it completes
    build = build or find_build(sha)
    if not build:
        return

    # The report is missing if MultiQC found no matching test data; the screenshot if it failed
    artifacts = {a.name: a for a in build.get_artifacts()}
    report = artifacts.get(f"multiqc_report_pr{pr.number}.html")
    shot = Path(f"report_pr{pr.number}.png")
    shot_artifact = artifacts.get(shot.name)
    shot.unlink(missing_ok=True)
    if shot_artifact:
        # archive: false artifacts download as the raw file
        resp = requests.get(
            shot_artifact.archive_download_url,
            headers={"Authorization": f"Bearer {os.environ['GH_TOKEN']}"},
            timeout=60,
        )
        if resp.ok:
            shot.write_bytes(resp.content)

    for c in mine:
        try:
            c.minimize()
        except Exception as e:
            print(f"::warning::Couldn't hide comment {c.id} on PR #{pr.number}: {e}")

    mods = ", ".join(f"`{m}`" for m in modules)
    run_link = f" See [run]({build.html_url}).<!-- preview-sha: {sha} -->"
    if build.conclusion != "success":
        post(pr, f"❌ {MARKER} failed.{run_link}")
    elif not report:
        post(pr, f"⚠️ No {MARKER} generated for {mods}.{run_link}")
    else:
        ok = f"✅ **[{MARKER}]({build.html_url}/artifacts/{report.id})** generated for {mods}"
        details = f"\n\n<details><summary>Screenshot</summary>\n\n![MultiQC report screenshot](./{shot})\n\n</details>"
        if not (shot_artifact and shot.exists()):
            post(pr, f"{ok} (⚠️ screenshot failed).{run_link}")
        else:
            try:
                post(pr, f"{ok}.{run_link}{details}", attach=shot)
            except subprocess.CalledProcessError:
                shot_url = f"{build.html_url}/artifacts/{shot_artifact.id}"
                post(pr, f"{ok} (⚠️ screenshot upload failed, [view it here]({shot_url})).{run_link}")
    react("rocket")


def main() -> None:
    react("eyes")
    prs, rule, build = get_targets()
    failed = False
    for pr in prs:
        try:
            process(pr, rule, build)
        except Exception as e:
            print(f"::warning::Preview comment for PR #{pr.number} failed: {e}")
            failed = True
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
