---
name: reviewing-prs
description: |
  Review MultiQC pull requests or local branch diffs — scan changes for anti-patterns against the project's rules, check PR hygiene (issue link, test data, screenshots), and produce a structured comment with severity-tagged findings. Use when the user asks to review a MultiQC PR, asks for feedback on a branch/diff, runs `gh pr view`/`gh pr diff`, or wants a code review of MultiQC changes (new modules, plot tweaks, core refactors, docs, anything).
---

# Review MultiQC Pull Requests

## Two modes

This skill runs in two contexts. Detect which one before deciding whether to
post:

- **CI / automated**: invoked by `.github/workflows/claude-code-review.yml`
  (triggers on PR opened or `/review` comment). The workflow runs in a
  GitHub Action, expects the agent to post the review via `gh pr comment`
  on its own, and provides a `PR_NUMBER` env var. In this mode: post the
  comment yourself, do not ask for confirmation.
- **Interactive**: the user is asking for a review locally (e.g. "review
  this branch", `gh pr view <n>` followed by "what do you think"). In
  this mode: produce the review text, show it to the user, and wait
  before posting. The user has a memory entry asking not to push or
  comment without explicit instruction.

Heuristic: if `$PR_NUMBER` is set in the environment, or the prompt
contains "MUST post" / "post the comment", you're in CI mode.

## Workflow

1. **Fetch the change**: `gh pr diff <n>`, `gh pr view <n>`, or
   `git diff main...HEAD` if reviewing a local branch. Pull the PR
   description too — half the review is checking the description, not
   just the code.
2. **Classify** the change: new module, plot/core change, bug fix, docs,
   refactor, dependency bump. The review focus shifts based on type.
3. **Scan the diff** for known anti-patterns (see Rule checklist below).
   Don't read files outside the diff unless an in-diff change clearly
   depends on context.
4. **Check PR hygiene** (see below).
5. **Write up** using the Output template.
6. **Deliver** based on mode: CI → `gh pr comment <n> --body-file -` with
   a heredoc. Interactive → print to chat and wait.

## Severity scheme

Tag each finding so the author knows what's required vs optional:

| Tag               | Meaning                                                                                                  |
| ----------------- | -------------------------------------------------------------------------------------------------------- |
| 🚫 **Blocker**    | Must fix before merge: breaks lint/tests/CI, contradicts a mandatory rule, or produces incorrect output. |
| 🔧 **Change**     | Should fix: clear deviation from project conventions or a real bug, but report would still render.       |
| 💡 **Suggestion** | Consider this: a better approach exists but the current code is reasonable.                              |
| 🪶 **Nit**        | Tiny / style preference. Skip if the author seems pressed for time.                                      |

Always cite the file and line: `multiqc/modules/foo/foo.py:123`.

## Rule checklist

The rules to check the diff against live in the existing skill files — do not
duplicate them here. Load whichever apply to the change type:

- **All MultiQC code** → [../../../CLAUDE.md](../../../CLAUDE.md) (style, em-dash
  ban, crash-loudly, helper-vs-inline, mandatory module rules).
- **New module / module change** → [../implementing-new-modules/SKILL.md](../implementing-new-modules/SKILL.md)
  "Common Pitfalls" section is the primary review checklist. Also load
  [../implementing-new-modules/implementation-checklist.md](../implementing-new-modules/implementation-checklist.md)
  for the `must`-marked calls and section-alert / colour-scale guidance.
- **Plot/template change** → check that any new IDs are unique (the `--lint` CI
  job catches HTML ID collisions).
- **Docs change** → confirm generic docs use placeholder names, not specific tool
  names. No em-dashes.

## PR hygiene

For a module PR:

- 🚫 `Closes #XXXX` in the PR description (links the originating issue).
- 🚫 Test data referenced — either added to the `MultiQC/test-data` repo or
  linked. Empty test runs aren't a review.
- 🔧 PR body has a brief summary, then a `<details>` block for the full
  write-up.
- 🔧 Sample report screenshot for any visible UI / plot change.
- 🔧 No em-dashes in PR title, body, or commit messages.

For a bug fix:

- 🔧 A regression test exists (or rationale for why it can't).
- 💡 Commit message says **why**, not just **what**.

For any PR:

- 🚫 CI is green (don't review around failing CI without flagging it).

## Things NOT to flag

Skip these to keep the review signal-to-noise high:

- Anything ruff, mypy, or `code_checks.py` already enforces — CI handles it.
- Pure stylistic taste with no rule behind it.
- Pre-existing issues outside the diff (unless the diff makes them worse).
- Speculative "what if a future tool added X" — react to what's there.
- The author's choice of variable names, unless genuinely confusing.

## Output template

Brief summary above the fold, full review wrapped in `<details>` — this
matches the CI workflow's required format and works fine for interactive
display too:

```markdown
<!-- 1–3 sentence summary of what the PR does and overall verdict -->

<details>
<summary>Full review</summary>

### Findings

🚫 **Blocker** — `path/to/file.py:42`

<!-- What's wrong, why, suggested fix -->

🔧 **Change** — `path/to/file.py:120`

<!-- ... -->

💡 **Suggestion** — `path/to/file.py:200`

<!-- ... -->

🪶 **Nit** — `path/to/file.py:301`

<!-- ... -->

### What's good

<!-- 2–3 bullets acknowledging strengths. Be specific. -->

</details>
```

Cap the review at ~10 findings unless the change is genuinely large;
otherwise prioritise blockers and changes. Author can ask for more.
