"""
Post MultiQC preview report comments on PRs, from builds made by module-report-build.yml.

Run by .github/workflows/module-report.yml. Never runs PR code: it only reads build
results and artifacts through the GitHub API, and comments with the bot token.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, List, Optional, Tuple

REPO = os.environ["REPO"]
EVENT = os.environ["EVENT"]
RUN_URL_BASE = f"{os.environ['GITHUB_SERVER_URL']}/{REPO}/actions/runs"
IS_REPORT = bool(re.search(r"(?m)^/report", os.environ.get("COMMENT_BODY", "")))

# Every bot comment contains this, so it can be found and hidden later
MARKER = "MultiQC preview report"
SHA_RE = re.compile(r"<!-- preview-sha: ([0-9a-f]+) -->")
# Keep in sync with module-report-build.yml
MODULE_RE = re.compile(r"^multiqc/modules/([A-Za-z0-9_]+)/", re.MULTILINE)

HOW = """

<details><summary>How this works</summary>

The bot posts a preview report on PRs that edit modules:

1. When the first build finishes, a few minutes after the PR opens. If you push more fixes before then, the report uses the latest commit.
2. Each night, if there are new commits since the last report.
3. When anyone comments `/report`. If the build is still running, the report posts when it finishes.

The bot only posts a new report when there are new commits. If there are none, it reacts 😕 to your comment. The bot hides older report comments.

</details>"""

COMMENTS_QUERY = """
query($owner: String!, $name: String!, $pr: Int!) {
  repository(owner: $owner, name: $name) { pullRequest(number: $pr) {
    comments(last: 100) { nodes { id isMinimized viewerDidAuthor body } } } } }"""

MINIMIZE_MUTATION = """
mutation($id: ID!) { minimizeComment(input: {subjectId: $id, classifier: OUTDATED}) { clientMutationId } }"""


def gh(*args: str, bot: bool = False) -> str:
    """Run gh. bot=True uses the bot PAT: needed for writes, viewerDidAuthor and --attach."""
    env = dict(os.environ, GH_TOKEN=os.environ["BOT_TOKEN"]) if bot else None
    return subprocess.run(["gh", *args], env=env, check=True, stdout=subprocess.PIPE, text=True).stdout


def gh_json(*args: str, bot: bool = False) -> Any:
    return json.loads(gh(*args, bot=bot))


def react(content: str) -> None:
    if IS_REPORT:
        gh(
            "api",
            f"repos/{REPO}/issues/comments/{os.environ['COMMENT_ID']}/reactions",
            "-f",
            f"content={content}",
            bot=True,
        )


def pr_head(pr: str) -> str:
    return gh("pr", "view", pr, "--repo", REPO, "--json", "headRefOid", "-q", ".headRefOid").strip()


def post(pr: str, body: str, attach: Optional[str] = None) -> None:
    args = ["pr", "comment", pr, "--repo", REPO, "--body", body + HOW]
    if attach:
        args += ["--attach", attach]
    gh(*args, bot=True)


def get_targets() -> Tuple[List[Tuple[str, str]], str, Optional[Tuple[str, str]]]:
    """
    Return (PR, head SHA) pairs, the posting rule, and the build (run ID, conclusion) if known.

    Rule "first": post only if the PR has no report yet, or has an unanswered /report.
    Rule "new": post only if the head moved since the last report.
    """
    if EVENT == "issue_comment":
        pr = os.environ["ISSUE"]
        return [(pr, pr_head(pr))], "new", None
    if EVENT == "schedule":
        prs = gh_json("pr", "list", "--repo", REPO, "--state", "open", "--limit", "500", "--json", "number,headRefOid")
        return [(str(p["number"]), p["headRefOid"]) for p in prs], "new", None

    # workflow_run. Fork runs have an empty pull_requests payload, so the build's
    # run-name carries the PR number.
    title = os.environ["WR_TITLE"]
    match = re.match(r"#(\d+) ", title)
    if not match:
        sys.exit(f"Can't parse PR from run name: {title}")
    pr = match.group(1)
    sha = pr_head(pr)
    build = (os.environ["WR_ID"], os.environ["WR_CONCLUSION"])
    if os.environ["WR_EVENT"] != "pull_request":
        return [(pr, sha)], "new", build
    # Also rejects a fork run naming someone else's PR
    if sha != os.environ["WR_SHA"]:
        print(f"PR #{pr} head moved on or doesn't match this build")
        return [], "first", build
    return [(pr, sha)], "first", build


def find_build(sha: str) -> Optional[Tuple[str, str]]:
    runs = gh_json(
        "api", f"repos/{REPO}/actions/workflows/module-report-build.yml/runs?head_sha={sha}&status=completed"
    )["workflow_runs"]
    for run in runs:
        if run["conclusion"] in ("success", "failure"):
            return str(run["id"]), run["conclusion"]
    return None


def process(pr: str, sha: str, rule: str, build: Optional[Tuple[str, str]]) -> None:
    modules = sorted(set(MODULE_RE.findall(gh("pr", "diff", pr, "--repo", REPO, "--name-only"))))
    if not modules:
        if IS_REPORT:
            post(pr, f"ℹ️ No edited modules found, so no {MARKER}.")
        return

    owner, name = REPO.split("/")
    nodes = gh_json(
        "api",
        "graphql",
        "-F",
        f"owner={owner}",
        "-F",
        f"name={name}",
        "-F",
        f"pr={pr}",
        "-f",
        f"query={COMMENTS_QUERY}",
        "--jq",
        ".data.repository.pullRequest.comments.nodes",
        bot=True,
    )
    mine = [i for i, n in enumerate(nodes) if n["viewerDidAuthor"] and MARKER in n["body"]]
    shas = [m.group(1) for i in mine if (m := SHA_RE.search(nodes[i]["body"]))]
    last_sha = shas[-1] if shas else ""
    # A /report that arrived while the build was still running, not yet answered by a report
    pending = any(re.search(r"(?m)^/report", n["body"]) for n in nodes[(mine[-1] + 1 if mine else 0) :])

    if last_sha == sha:
        react("confused")
        return
    if rule == "first" and last_sha and not pending:
        return

    # No finished build yet: its workflow_run posts the report when it completes
    build = build or find_build(sha)
    if not build:
        return
    run_id, result = build
    run_url = f"{RUN_URL_BASE}/{run_id}"

    # The report is missing if MultiQC found no matching test data; the screenshot if it failed
    artifacts = {
        a["name"]: a["id"] for a in gh_json("api", f"repos/{REPO}/actions/runs/{run_id}/artifacts")["artifacts"]
    }
    report_id = artifacts.get(f"multiqc_report_pr{pr}.html")
    shot = Path(f"report_pr{pr}.png")
    shot_id = artifacts.get(shot.name)
    shot.unlink(missing_ok=True)
    if shot_id:
        # archive: false artifacts download as the raw file
        with shot.open("wb") as fh:
            if subprocess.run(["gh", "api", f"repos/{REPO}/actions/artifacts/{shot_id}/zip"], stdout=fh).returncode:
                shot.unlink()

    for i in mine:
        if not nodes[i]["isMinimized"]:
            try:
                gh("api", "graphql", "-f", f"id={nodes[i]['id']}", "-f", f"query={MINIMIZE_MUTATION}", bot=True)
            except subprocess.CalledProcessError:
                print(f"::warning::Couldn't hide comment {nodes[i]['id']} on PR #{pr}")

    mods = ", ".join(f"`{m}`" for m in modules)
    ok = f"✅ **[{MARKER}]({run_url}/artifacts/{report_id})** generated for {mods}"
    run_link = f" See [run]({run_url}).<!-- preview-sha: {sha} -->"
    if result != "success":
        body = f"❌ {MARKER} failed."
    elif not report_id:
        body = f"⚠️ No {MARKER} generated for {mods}."
    elif shot.exists():
        details = f"\n\n<details><summary>Screenshot</summary>\n\n![MultiQC report screenshot](./{shot})\n\n</details>"
        try:
            post(pr, f"{ok}.{run_link}{details}", attach=f"./{shot}")
            react("rocket")
            return
        except subprocess.CalledProcessError:
            body = f"{ok} (⚠️ screenshot upload failed, [view it here]({run_url}/artifacts/{shot_id}))."
    else:
        body = f"{ok} (⚠️ screenshot failed)."
    post(pr, body + run_link)
    react("rocket")


def main() -> None:
    react("eyes")
    targets, rule, build = get_targets()
    failed = False
    for pr, sha in targets:
        try:
            process(pr, sha, rule, build)
        except Exception as e:
            print(f"::warning::Preview comment for PR #{pr} failed: {e}")
            failed = True
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
